#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NCBI_DOWNLOAD } from '../modules/ncbi-dl.nf'

workflow GENOME_DOWNLOAD {
    main:
    // Download the files
    NCBI_DOWNLOAD(params.genome)
    ReferenceFasta = NCBI_DOWNLOAD.out.fasta
    ReferenceGenbank = NCBI_DOWNLOAD.out.genbank

    // Measure the size of the genome based on the genbank record
    GenomeSize = ReferenceGenbank
        .splitText()
        .first()
        .map{ l -> l.replaceAll(/\s+/, ',') }
        .splitCsv()
        .flatten()
        .take(3)
        .last()

    // Get the name of the reference genome
    ReferenceName = ReferenceGenbank
        .splitText()
        .take(2)
        .last()
        .map{ s -> s.replace('DEFINITION  ', '') }
        .map{ s -> s.replace(',', '') }
        .map{ s -> s.replace('.', '') }
        .map{ s -> s.replace(' ', '_') }
        .map{ s -> s.replace('/', '_') }
        .map{ s -> s.trim() }

    // Process the files
    INDEXING(ReferenceFasta, ReferenceName)
    IndexedReference = INDEXING.out
    ANNOTATION(ReferenceGenbank, ReferenceName)
    AnnotatedReference = ANNOTATION.out

    emit:
    indexedFasta = IndexedReference
    referenceAnnotations = AnnotatedReference
    genomesize = GenomeSize
}

// Index the reference genome for use with Samtools
process INDEXING {
    label 'samtools'
    label 'process_low'
    publishDir "${params.outdir}/reference", mode: "${params.publish_dir_mode}"

    input:
    file(reference)
    val(refname)

    output:
    tuple file("${refname}.fasta"), file("*.fai")

    script:
    """
    # Create a reference genome index
    cp ${reference} ${refname}.fasta
    samtools faidx ${refname}.fasta
    """
}

// Process the reference genome's feature table into GFF format
process ANNOTATION {
    label 'seqret'
    label 'process_low'
        publishDir "${params.outdir}/reference", mode: "${params.publish_dir_mode}"

    input:
    file(reference)
    val(refname)

    output:
    file "${refname}.gff"

    shell:
    '''
    seqret -sequence !{reference} -sformat1 genbank -feature -outseq ref.gff -osformat gff -auto
    head -n $(($(grep -n '##FASTA' ref.gff | cut -d : -f 1) - 1)) ref.gff > !{refname}.gff
    '''
}
