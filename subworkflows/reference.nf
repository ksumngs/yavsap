#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EFETCH } from '../modules/local/modules/efetch/main.nf'
include { ESEARCH } from '../modules/local/modules/esearch/main.nf'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/modules/samtools/faidx/main.nf'

workflow GENOME_DOWNLOAD {
    main:
    Channel
        .of([
                ['id': 'reference', 'single_end': null, 'strandedness': null],
                'nucleotide',
                params.genome
            ])
        .set{ SearchParameters }

    ESEARCH(SearchParameters)

    ESEARCH.out.xml
        .combine(
            Channel.of(
                ['fasta', null]
            )
        )
        .set{ FetchParameters }

    EFETCH(FetchParameters)

    SAMTOOLS_FAIDX(
        EFETCH.out.txt
    )

    EFETCH.out.txt
        .join(SAMTOOLS_FAIDX.out.fai)
        .set{ fasta }

    emit:
    fasta
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
