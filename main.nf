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

include { reference_genome_pull } from './subworkflows/reference.nf'
include { trimming }              from './subworkflows/trimming.nf'
include { assembly }              from './subworkflows/assembly.nf'
include { read_filtering }        from './subworkflows/filtering.nf'
include { haplotyping }           from './subworkflows/haplotype.nf'
include { simulated_reads }       from './test/test.nf'

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
    reference_genome_pull()
    IndexedReference = reference_genome_pull.out.indexedreference
    AnnotatedReference = reference_genome_pull.out.annotatedreference
    GenomeSize = reference_genome_pull.out.genomesize

    // Bring in the reads files
    if (params.sra) {
        simulated_reads()
        RawReads = simulated_reads.out
    }
    else {
        if (params.ont) {
            RawReads = Channel
                .fromPath("${params.input}/*.{fastq,fq}.gz")
                .map{ file -> tuple(file.simpleName, file) }
        }
        else {
            RawReads = Channel
                .fromFilePairs("${params.input}/*{R1,R2,_1,_2}*.{fastq,fq}.gz")
        }
    }

    RawReads | sample_rename | trimming

    if (params.ont) {
        sample_rename.out | nanostat
        QcReports = nanostat.out
    }
    else {
        sample_rename.out | interleave | fastqc
        QcReports = fastqc.out
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
        haplotyping(Alignments, Assemblies, IndexedReference, AnnotatedReference)
        PhyloTrees = haplotyping.out
    }
    else {
        PhyloTrees = Channel.from([])
    }

    multiqc(KrakenReports
        .concat(trimming.out.report)
        .concat(QcReports)
        .collect())

    // Put a pretty bow on everything
    presentation_generator()
}

process sample_rename {
    label 'process_low'

    input:
    tuple val(givenName), file(readsFiles)

    output:
    tuple val(sampleName), file("out/*.fastq.gz")

    script:
    sampleName = givenName.split('_')[0]
    if (params.ont) {
        """
        mkdir out
        mv ${readsFiles[0]} out/${sampleName}.fastq.gz
        """
    }
    else {
        """
        mkdir out
        mv ${readsFiles[0]} out/${sampleName}_R1.fastq.gz
        mv ${readsFiles[1]} out/${sampleName}_R2.fastq.gz
        """
    }
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
    publishDir "${params.outdir}/alignment", mode: "${params.publish_dir_mode}"

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
    file '*'

    output:
    path 'multiqc_report.html' optional true
    path 'multiqc_data'        optional true

    script:
    """
    cp ${workflow.projectDir}/multiqc_config.yaml .
    multiqc .
    """
}

// Create a viewer of all the assembly files
process presentation_generator {
    label 'process_low'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    output:

    file '_css/*.css'
    file '_views/*.pug'
    file 'index.js'
    file 'package.json'
    file 'package-lock.json'
    file 'favicon.ico'

    script:
    """
    cp -r ${workflow.projectDir}/visualizer/{_css,_views,index.js,package.json,package-lock.json,favicon.ico} .
    """
}
