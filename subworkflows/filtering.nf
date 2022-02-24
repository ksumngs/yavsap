#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KRAKEN2 } from '../modules/ksumngs/nf-modules/kraken2/main.nf'
include { KRAKENTOOLS_EXTRACT } from '../modules/ksumngs/nf-modules/krakentools/extract/main.nf'

workflow FILTERING {
    take:
    reads
    kraken2_db
    filter

    main:
    KRAKEN2(reads, kraken2_db)

    KRAKEN2.out.kreport.set{ log_out }

    if (filter == 'classified') {
        KRAKEN2.out.classified.set{ filtered }
    }
    else if ( filter == 'unclassified') {
        KRAKEN2.out.unclassified.set{ filtered }
    }
    else {
        KRAKENTOOLS_EXTRACT(
            reads
                .join(KRAKEN2.out.kraken)
                .join(KRAKEN2.out.kreport),
            filter
        )
        KRAKENTOOLS_EXTRACT.out.fastq.set{ filtered }
    }

    emit:
    filtered
    log_out
}
