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
include { CLOSEST_REFERENCE } from '../subworkflows/local/closest-reference'
include { CONSENSUS } from '../subworkflows/local/consensus'
include { FILTERING } from '../subworkflows/local/filtering.nf'
include { GENOME_DOWNLOAD } from '../subworkflows/local/genomes'
include { HAPLOTYPING } from '../subworkflows/local/haplotype.nf'
include { PHYLOGENETIC_TREE } from '../subworkflows/local/phylogenetics.nf'
include { QC } from '../subworkflows/local/qc.nf'
include { READS_INGEST } from '../subworkflows/local/ingest.nf'
include { REFERENCE_DOWNLOAD } from '../subworkflows/local/reference'
include { TRIMMING } from '../subworkflows/local/trimming.nf'
include { VARIANTS } from '../subworkflows/local/variants'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main.nf'
include { KRAKEN2_DBPREPARATION } from '../modules/local/kraken2/dbpreparation.nf'
include { MINIMAP2_ALIGN as MINIMAP2_REALIGN } from '../modules/ksumngs/nf-modules/minimap2/align/main'
include { MINIMAP2_ALIGN } from '../modules/nf-core/modules/minimap2/align/main'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main.nf'
include { PRESENTATION } from '../subworkflows/presentation.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_REINDEX } from '../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main'

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

    //
    // SUBWORKFLOW: Download reference genome from NCBI
    //
    REFERENCE_DOWNLOAD()
    REFERENCE_DOWNLOAD.out.fasta.set{ ch_reference_fasta }
    ch_versions = ch_versions.mix(REFERENCE_DOWNLOAD.out.versions)

    //
    // MODULE: Align reads into BAM format using minimap2
    //
    MINIMAP2_ALIGN(ch_filtered, ch_reference_fasta, true, false, false)
    MINIMAP2_ALIGN.out.bam.set{ ch_bam }
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    //
    // MODULE: Index BAM reads using Samtools
    //
    SAMTOOLS_INDEX(ch_bam)
    SAMTOOLS_INDEX.out.bai.set{ ch_bai }
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // SUBWORKFLOW: Find the consensus sequence
    //
    CONSENSUS(ch_bam, ch_bai, ch_reference_fasta)
    CONSENSUS.out.fasta.set{ ch_consensus_fasta }
    ch_versions = ch_versions.mix(CONSENSUS.out.versions)

    //
    // SUBWORKFLOW: Download and reformat the strain reference genomes
    //
    GENOME_DOWNLOAD("${params.genome_list}")
    GENOME_DOWNLOAD.out.strain.set{ ch_genome_strain }
    GENOME_DOWNLOAD.out.fasta.set{ ch_genome_fasta }
    ch_versions = ch_versions.mix(GENOME_DOWNLOAD.out.versions)

    //
    // SUBWORKFLOW: Find the closest strain to each consensus sequence
    //
    CLOSEST_REFERENCE(ch_consensus_fasta, ch_genome_strain, ch_genome_fasta)
    CLOSEST_REFERENCE.out.accession.set{ ch_closest_accession }
    CLOSEST_REFERENCE.out.strain.set{ ch_closest_strain }
    CLOSEST_REFERENCE.out.fasta.set{ ch_closest_reference }
    ch_versions = ch_versions.mix(CLOSEST_REFERENCE.out.versions)

    //
    // MODULE: Realign reads into BAM format using minimap2
    //
    MINIMAP2_REALIGN(ch_filtered.join(ch_closest_reference), true, false, false)
    MINIMAP2_REALIGN.out.bam.set{ ch_realigned_bam }
    ch_versions = ch_versions.mix(MINIMAP2_REALIGN.out.versions.first())

    //
    // MODULE: Index new BAM reads using Samtools
    //
    SAMTOOLS_REINDEX(ch_realigned_bam)
    SAMTOOLS_REINDEX.out.bai.set{ ch_realigned_bai }
    ch_versions = ch_versions.mix(SAMTOOLS_REINDEX.out.versions.first())

    //
    // SUBWORKFLOW: Variant calling
    //
    VARIANTS(ch_realigned_bam, ch_realigned_bai, ch_closest_reference)
    VARIANTS.out.vcf.set{ ch_vcf }
    ch_versions = ch_versions.mix(VARIANTS.out.versions.first())

    //
    // SUBWORKFLOW: Haplotype calling
    //
    ch_haplotype_fasta = Channel.empty()
    ch_haplotype_yaml = Channel.empty()
    if (!params.skip_haplotype) {
        HAPLOTYPING(ch_realigned_bam, ch_vcf, ch_closest_reference)
        HAPLOTYPING.out.fasta.set{ ch_haplotype_fasta }
        HAPLOTYPING.out.yaml.set { ch_haplotype_yaml }
        ch_versions = ch_versions.mix(HAPLOTYPING.out.versions)
    }

    //
    // SUBWORKFLOW: Phylogenetics
    //
    PHYLOGENETIC_TREE(
        ch_haplotype_fasta,
        ch_consensus_fasta,
        ch_genome_fasta,
        ch_genome_strain
    )
    PHYLOGENETIC_TREE.out.tree.set{ ch_tree }
    ch_versions = ch_versions.mix(PHYLOGENETIC_TREE.out.versions)

    //
    // MODULE: Get the versions of each bioinformatics tool
    //
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

    VersionFiles = VersionFiles.mix(CLOSEST_REFERENCE.out.versions)

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
