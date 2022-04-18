#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ALIGNMENT } from './subworkflows/alignment.nf'
include { CLOSEST_REFERENCE } from './subworkflows/closest-reference.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/modules/custom/dumpsoftwareversions/main.nf'
include { FILTERING } from './subworkflows/filtering.nf'
include { GENOME_DOWNLOAD } from './subworkflows/reference.nf'
include { HAPLOTYPING } from './subworkflows/haplotype.nf'
include { KRAKEN2_DBPREPARATION } from './modules/local/kraken2/dbpreparation.nf'
include { MULTIQC } from './modules/nf-core/modules/multiqc/main.nf'
include { PHYLOGENETIC_TREE } from './subworkflows/phylogenetics.nf'
include { PRESENTATION } from './subworkflows/presentation.nf'
include { QC } from './subworkflows/qc.nf'
include { READS_INGEST } from './subworkflows/ingest.nf'
include { TRIMMING } from './subworkflows/trimming.nf'
include { cowsay } from './lib/cowsay.nf'
include { yavsap_logo } from './lib/logo.nf'

cowsay(yavsap_logo())

if (params.help) {
    log.info(
        """\
        YAVSAP (Yet Another Viral Subspecies Analysis Pipeline) - Intra-sample viral
        population analysis

        Usage:

            nextflow run ksumngs/yavsap

        Options:

            --input             Relative or absolute path to directory containing
                                gzipped fastq files
                                    type: path, default: .

            --platform          Type of reads to process. Options are 'illumina' and
                                'nanopore'
                                    type: string, default: none

            --genome            NCBI accession number of the reference genome to
                                align reads against
                                    type: string, default: 'NC_001437.1'

            --kraken2_db        Kraken2-compatible database for classifying reads
                                    type: path, default: none

            --keep_taxid        Space-separated list of NCBI taxids to keep and
                                process after classifying
                                    type: string, default: '0 10239'

            --outdir            Directory in which to place results
                                    type: path, default: ./results

            --help              Print this message and exit

        For more information on usage and parameters, visit the website at
            https://ksumngs.github.io/yavsap
    """.stripIndent()
    )

    exit 0
}

if (params.platform != 'illumina' && params.platform != 'nanopore') {
    log.error "ERROR: --platform <illumina,nanopore> must be specified"
    exit 1
}

log.info(
    """\
    Input folder:           ${params.input}
    Sequencing platform:    ${params.platform}
    Reference genome:       ${params.genome}
    Kraken2 Database:       ${params.kraken2_db}
    Taxonomic Ids:          '${params.keep_taxid}'
    Output folder           ${params.outdir}
    Diagnostics folder:     ${params.tracedir}
    """.stripIndent()
)

workflow {
    LogFiles = Channel.empty()
    VersionFiles = Channel.empty()

    GENOME_DOWNLOAD()
    ReferenceGenome = GENOME_DOWNLOAD.out.fasta
    VersionFiles = VersionFiles.mix(GENOME_DOWNLOAD.out.versions)

    // Bring in the reads files
    READS_INGEST()
    RawReads = READS_INGEST.out.sample_info
    VersionFiles = VersionFiles.mix(READS_INGEST.out.versions)

    if (!params.skip_qc) {
        QC(RawReads)
        LogFiles = LogFiles.mix(QC.out.report)
        VersionFiles = VersionFiles.mix(QC.out.versions)
    }

    if (!params.skip_trimming) {
        TRIMMING(RawReads)
        TRIMMING.out.fastq.set{ TrimmedReads }
        LogFiles = LogFiles.mix(TRIMMING.out.log_out)
        VersionFiles = VersionFiles.mix(TRIMMING.out.versions)
    }
    else {
        RawReads.set { TrimmedReads }
    }

    KronaChart = Channel.of([])

    if (!params.skip_filtering) {
        KrakenDb = file("${params.kraken2_db}", checkIfExists: true)
        if (KrakenDb.isDirectory()) {
            // This is a local database, and everything is ready to pass to the
            // filtering process

        }
        else {
            if (KrakenDb.getExtension() == 'k2d') {
                // The user got confused, and passed a database file, we'll try to
                // correct it for them
                log.warn "WARNING: ${params.kraken2_db} appears to be a file that is a *part* of a Kraken2 database."
                log.warn "         Kraken databases are folders that contain multiple files."
                log.warn "         YAVSAP will attempt to use the parent directory as the database, but it might fail!"
                KrakenDb = KrakenDb.getParent()
            }
            else {
                // We'll assume this is a tarballed database
                KRAKEN2_DBPREPARATION(KrakenDb)
                KrakenDb = KRAKEN2_DBPREPARATION.out.db
                VersionFiles = VersionFiles.mix(KRAKEN2_DBPREPARATION.out.versions)
            }
        }
        FILTERING(
            TrimmedReads,
            KrakenDb,
            "${params.keep_taxid}"
        )
        FILTERING.out.filtered.set { FilteredReads }
        FILTERING.out.krona.set { KronaChart }
        LogFiles = LogFiles.mix(FILTERING.out.log_out)
        VersionFiles = VersionFiles.mix(FILTERING.out.versions)
    }
    else {
        TrimmedReads.set { FilteredReads }
    }

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
        ReferenceGenome,
        genomeFile
    )

    VersionFiles = VersionFiles.mix(CLOSEST_REFERENCE.out.versions)

    if (!params.skip_haplotype) {
        HAPLOTYPING(
            CLOSEST_REFERENCE.out.bam
                .join(
                    CLOSEST_REFERENCE.out.bai
                ),
            CLOSEST_REFERENCE.out.fasta
        )

        VersionFiles = VersionFiles.mix(HAPLOTYPING.out.versions)

        if (!params.skip_phylogenetics) {
            PHYLOGENETIC_TREE(
                HAPLOTYPING.out.fasta,
                CLOSEST_REFERENCE.out.consensus_fasta,
                CLOSEST_REFERENCE.out.genome_fasta,
                genomeFile
            )

            VersionFiles = VersionFiles.mix(PHYLOGENETIC_TREE.out.versions)
        }
    }

    LogFiles = LogFiles
        .map{ (it instanceof Path) ? it : it.drop(1) }
        .mix(Channel.of(file("${workflow.projectDir}/assets/multiqc_config.yml")))
        .flatten()
        .collect()

    MULTIQC(LogFiles)

    PRESENTATION(
        ALIGNMENT.out.bam,
        ALIGNMENT.out.bai,
        GENOME_DOWNLOAD.out.fasta,
        GENOME_DOWNLOAD.out.fai,
        CLOSEST_REFERENCE.out.consensus_fasta,
        CLOSEST_REFERENCE.out.accession,
        CLOSEST_REFERENCE.out.strain,
        HAPLOTYPING.out.yaml,
        HAPLOTYPING.out.fasta
    )

    CUSTOM_DUMPSOFTWAREVERSIONS(VersionFiles.unique().collectFile(name: 'collated_versions.yml'))

    // Put a pretty bow on everything
    PRESENTATION_GENERATOR()
}

// Create a viewer of all the assembly files
process PRESENTATION_GENERATOR {
    label 'process_low'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    output:
    file '_css/*.css'
    file '_js/*.js'
    file '_views/*.pug'
    file 'index.js'
    file 'package.json'
    file 'package-lock.json'
    file 'favicon.ico'

    script:
    """
    cp -r ${workflow.projectDir}/visualizer/{_css,_js,_views,index.js,package.json,package-lock.json,favicon.ico} .
    """
}
