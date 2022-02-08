#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { cowsay } from './lib/cowsay.nf'

if (params.help) {
    cowsay(
    """\
====================================================================================
                                     YAVSAP
                (Yet Another Viral Subspecies Analysis Pipeline)
====================================================================================

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

        --genome            NCBI accession number of the reference genome to align
                            reads against
                                type: string, default: 'NC_001437.1'

        --kraken2_db        Kraken2-compatible database for classifying reads
                                type: path, default: none

        --keep_taxid        Space-separated list of NCBI taxids to keep and process
                            after classifying
                                type: string, default: '0 10239'

        --outdir            Directory in which to place results
                                type: path, default: ./results

        --help              Print this message and exit

    For more information on usage and parameters, visit the website at
        https://ksumngs.github.io/yavsap
"""
)
exit 0
}

if (!params.ont && !params.pe) {
    log.error "ERROR: --platform <illumina,nanopore> must be specified"
    exit 1
}

include { GENOME_DOWNLOAD } from './subworkflows/reference.nf'
include { READS_INGEST } from './modules/ingest.nf'
include { trimming }              from './subworkflows/trimming.nf'
include { assembly }              from './subworkflows/assembly.nf'
include { read_filtering }        from './subworkflows/filtering.nf'
include { haplotyping }           from './subworkflows/haplotype.nf'
include { SIMULATED_READS }       from './test/test.nf'

cowsay(
"""\
====================================================================================
                                     YAVSAP
                (Yet Another Viral Subspecies Analysis Pipeline)
====================================================================================

Input folder:           ${params.input}
Sequencing platform:    ${params.platform}
    Illumina?:          ${params.pe}
    Nanopore?:          ${params.ont}
Reference genome:       ${params.genome}
Kraken2 Database:       ${params.kraken2_db}
Taxonomic Ids:          '${params.keep_taxid}'
Output folder           ${params.outdir}
Diagnostics folder:     ${params.tracedir}
"""
)

workflow {
    GENOME_DOWNLOAD()
    IndexedReference = GENOME_DOWNLOAD.out.indexedFasta
    AnnotatedReference = GENOME_DOWNLOAD.out.referenceAnnotations
    GenomeSize = GENOME_DOWNLOAD.out.genomeSize

    // Bring in the reads files
    if (params.sra) {
        SIMULATED_READS()
        RawReads = SIMULATED_READS.out
    }
    else {
        READS_INGEST()
        RawReads = READS_INGEST.out
    }

    trimming(RawReads)

    if (params.skip_qc) {
        QcReport = Channel.from([])
    }
    else {
        if (params.ont) {
            InterleavedReads = RawReads
        }
        else {
            interleave(RawReads)
            InterleavedReads = interleave.out
        }

        fastqc(InterleavedReads)
        QcReport = fastqc.out
    }


    read_filtering(trimming.out.trimmedreads)
    KrakenReports = read_filtering.out.KrakenReports
    FilteredReads = read_filtering.out.FilteredReads

    if (!params.skip_assembly) {
        // _de novo_ assemble the viral reads
        assembly(FilteredReads, IndexedReference, GenomeSize)
        Assemblies = assembly.out.Contigs
        AlignedContigs = assembly.out.AlignedContigs
    }
    else {
        Assemblies = Channel.from([])
        AlignedContigs = Channel.from([])
    }

    // Realign reads to the reference genome
    reads_realign_to_reference(FilteredReads, IndexedReference)
    Alignments = reads_realign_to_reference.out

    AllAlignments = Alignments.join(AlignedContigs, remainder: true).flatMap{ n -> [n[1], n[2]] }.collect()

    if (!params.skip_haplotype) {
        haplotyping(FilteredReads, Alignments, IndexedReference, AnnotatedReference)
        //PhyloTrees = haplotyping.out
    }
    else {
        //PhyloTrees = Channel.from([])
    }

    MultiQCConfig = file("${workflow.projectDir}/multiqc_config.yaml")

    multiqc(MultiQCConfig,
        KrakenReports
        .concat(trimming.out.report)
        .concat(QcReport)
        .collect())

    // Put a pretty bow on everything
    presentation_generator()
}

process interleave {
    label 'seqtk'
    label 'process_low'

    input:
    tuple val(sampleName), path(readsFiles)

    output:
    tuple val(sampleName), path("${sampleName}.fastq.gz")

    script:
    """
    seqtk mergepe ${readsFiles} | gzip > ${sampleName}.fastq.gz
    """
}

process fastqc {
    label 'fastqc'
    label 'process_medium'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    path("${sampleName}_fastqc.zip")

    script:
    """
    fastqc -t ${task.cpus} ${readsFiles}
    """
}

process nanostat {
    label 'nanostat'
    label 'process_medium'

    input:
    tuple val(sampleName), file(readsFile)

    output:
    path("${sampleName}_nanostat.log")

    script:
    """
    NanoStat -t ${task.cpus} --fastq ${readsFile} > ${sampleName}_nanostat.log
    """
}

process reads_realign_to_reference {
    label 'minimap'

    input:
    tuple val(sampleName), file(readsFile)
    file(reference)

    output:
    tuple val(sampleName), file("${sampleName}.{bam,bam.bai}")

    script:
    minimapMethod = (params.pe) ? 'sr' : 'map-ont'
    """
    minimap2 -ax ${minimapMethod} -t ${task.cpus} --MD ${reference[0]} ${readsFile} | \
        samtools sort > ${sampleName}.bam
    samtools index ${sampleName}.bam
    """
}

process multiqc {
    label 'process_low'
    label 'multiqc'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    input:
    file(configFile)
    file '*'

    output:
    path 'multiqc_report.html' optional true
    path 'multiqc_data'        optional true

    script:
    """
    multiqc .
    """
}

// Create a viewer of all the assembly files
process presentation_generator {
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
