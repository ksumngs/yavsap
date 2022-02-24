#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NANOFILT } from '../modules/ksumngs/nf-modules/nanofilt/main.nf'
include { TRIMMOMATIC } from '../modules/ksumngs/nf-modules/trimmomatic/main.nf'

workflow TRIMMING {
    take:
    reads

    main:
    if (params.platform == 'illumina') {
        TRIMMOMATIC(reads)
        TRIMMOMATIC.out.fastq.set{ fastq }
        TRIMMOMATIC.out.log.set{ log_out }
    }
    else if (params.platform == 'nanopore') {
        NANOFILT(reads)
        NANOFILT.out.fastq.set{ fastq }
        NANOFILT.out.log.set{ log_out }
    }

    emit:
    fastq
    log_out
}
