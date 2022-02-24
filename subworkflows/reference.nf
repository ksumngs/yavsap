#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EDIRECT_EFETCH } from '../modules/ksumngs/nf-modules/edirect/efetch/main.nf'
include { EDIRECT_ESEARCH } from '../modules/ksumngs/nf-modules/edirect/esearch/main.nf'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/modules/samtools/faidx/main.nf'

workflow GENOME_DOWNLOAD {
    main:
    EDIRECT_ESEARCH("${params.genome}", 'nucleotide')
    EDIRECT_EFETCH(EDIRECT_ESEARCH.out.xml, 'fasta', '')
    EDIRECT_EFETCH.out.txt.set{ fasta }

    SAMTOOLS_FAIDX(
        fasta.map{ [
            ['id': 'reference'],
            it
        ] }
    )
    SAMTOOLS_FAIDX.out.fai.set{ fai }

    emit:
    fasta
    fai
}
