#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC } from '../modules/nf-core/modules/fastqc/main.nf'
include { NANOSTAT } from '../modules/ksumngs/nf-modules/nanostat/main.nf'

/// summary: |
///   Perform context-sensitive QC on fastq reads
workflow QC {
    take:
    reads

    main:
    if (params.platform == 'illumina') {
        FASTQC(reads)
        FASTQC.out.zip.set{ report }
    }
    else if (params.platform == 'nanopore') {
        NANOSTAT(reads)
        NANOSTAT.out.log.set{ report }
    }

    emit:
    report
}
