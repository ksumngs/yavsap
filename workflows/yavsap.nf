/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowYavsap.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2_db ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input does not exist!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FILTERING } from '../subworkflows/local/filtering.nf'
include { QC } from '../subworkflows/local/qc.nf'
include { READS_INGEST } from '../subworkflows/local/ingest.nf'
include { TRIMMING } from '../subworkflows/local/trimming.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { ALIGNMENT } from '../subworkflows/alignment.nf'
include { CLOSEST_REFERENCE } from '../subworkflows/closest-reference.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main.nf'
include { HAPLOTYPING } from '../subworkflows/haplotype.nf'
include { KRAKEN2_DBPREPARATION } from '../modules/local/kraken2/dbpreparation.nf'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main.nf'
include { PHYLOGENETIC_TREE } from '../subworkflows/phylogenetics.nf'
include { PRESENTATION } from '../subworkflows/presentation.nf'
include { REFERENCE_DOWNLOAD } from '../subworkflows/reference.nf'
include { cowsay } from '../lib/cowsay.nf'
include { yavsap_logo } from '../lib/logo.nf'

// include { MULTIQC } from '../modules/nf-core/modules/multiqc/main'
// include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow YAVSAP {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    READS_INGEST()
    READS_INGEST.out.sample_info.set{ ch_reads }
    ch_versions = ch_versions.mix(READS_INGEST.out.versions)

    //
    // SUBWORKFLOW: Run read QC
    //
    ch_qc = Channel.empty()
    if (!params.skip_qc) {
        QC(ch_reads)
        QC.out.report.set{ ch_qc }
        ch_versions = ch_versions.mix(QC.out.versions)
    }

    //
    // SUBWORKFLOW: Trim reads
    //
    ch_reads.set{ ch_trimmed }
    ch_trimlog = Channel.empty()
    if (!params.skip_trimming) {
        TRIMMING(ch_reads)
        TRIMMING.out.fastq.set{ ch_trimmed }
        TRIMMING.out.log_out.set{ ch_trimlog }
        ch_versions = ch_versions.mix(TRIMMING.out.versions)
    }

    //
    // SUBWORKFLOW: Kraken2 host read filtering
    //
    ch_trimmed.set{ ch_filtered }
    ch_krona = Channel.empty()
    ch_kreport = Channel.empty()
    if (!params.skip_filtering) {
        FILTERING(ch_trimmed, "${params.kraken2_db}", "${params.keep_taxid}")
        FILTERING.out.filtered.set{ ch_filtered }
        FILTERING.out.krona.set{ ch_krona }
        FILTERING.out.log_out.set{ ch_kreport }
        ch_versions = ch_versions.mix(FILTERING.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowYavsap.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

        LogFiles = Channel.empty()
    VersionFiles = Channel.empty()

    GENOME_DOWNLOAD()
    ReferenceGenome = GENOME_DOWNLOAD.out.fasta
    VersionFiles = VersionFiles.mix(GENOME_DOWNLOAD.out.versions)

    ALIGNMENT(FilteredReads, ReferenceGenome)
    ALIGNMENT.out.bam
        .join(ALIGNMENT.out.bai)
        .set { AlignedReads }

    VersionFiles = VersionFiles.mix(ALIGNMENT.out.versions)

    // Find the strain genomes list
    genomePath = params.genome_list
    genomeFile = file(genomePath, type: 'file')
    if (genomeFile.toFile().exists()) {
        genomeFile = [genomeFile]
    }
    else {
        genomePath = "${workflow.projectDir}/genomes/${params.genome_list}*"
        genomeFile = file(genomePath, checkIfExists: true, type: 'file')
    }

    // Realign reads to their closest strain
    CLOSEST_REFERENCE(
        ALIGNMENT.out.bam,
        ALIGNMENT.out.bai,
        ReferenceGenome,
        genomeFile
    )

    VersionFiles = VersionFiles.mix(CLOSEST_REFERENCE.out.versions)

    HaplotypeFastas = Channel.empty()
    HaplotypeYamls = Channel.empty()

    if (!params.skip_haplotype) {
        HAPLOTYPING(
            CLOSEST_REFERENCE.out.bam
                .join(
                    CLOSEST_REFERENCE.out.bai
                ),
            CLOSEST_REFERENCE.out.fasta
        )

        HAPLOTYPING.out.fasta.set{ HaplotypeFastas }
        HAPLOTYPING.out.yaml.set{ HaplotypeYamls }

        VersionFiles = VersionFiles.mix(HAPLOTYPING.out.versions)
    }

    PhylogeneticTree = Channel.of([])

    if (!params.skip_phylogenetics) {
        PHYLOGENETIC_TREE(
            HaplotypeFastas,
            CLOSEST_REFERENCE.out.consensus_fasta,
            CLOSEST_REFERENCE.out.genome_fasta,
            genomeFile
        )

        PHYLOGENETIC_TREE.out.tree.set{ PhylogeneticTree }

        VersionFiles = VersionFiles.mix(PHYLOGENETIC_TREE.out.versions)
    }

    LogFiles = LogFiles
        .map{ (it instanceof Path) ? it : it.drop(1) }
        .mix(Channel.of(file("${workflow.projectDir}/assets/multiqc_config.yml")))
        .flatten()
        .collect()

    // MULTIQC(LogFiles)
    VersionFiles = VersionFiles.mix(MULTIQC.out.versions)

    // Note: The Visualizer cannot be output if haplotyping is skipped
    PRESENTATION(
        ALIGNMENT.out.bam,
        ALIGNMENT.out.bai,
        GENOME_DOWNLOAD.out.fasta,
        GENOME_DOWNLOAD.out.fai,
        CLOSEST_REFERENCE.out.consensus_fasta,
        CLOSEST_REFERENCE.out.accession,
        CLOSEST_REFERENCE.out.strain,
        HaplotypeYamls,
        HaplotypeFastas,
        PhylogeneticTree,
        MULTIQC.out.report,
        KronaChart
    )
    VersionFiles = VersionFiles.mix(PRESENTATION.out.versions)

    // CUSTOM_DUMPSOFTWAREVERSIONS(VersionFiles.unique().collectFile(name: 'collated_versions.yml'))

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
