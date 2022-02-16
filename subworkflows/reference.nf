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
