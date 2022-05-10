#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KRAKEN2 } from '../modules/ksumngs/nf-modules/kraken2/main.nf'
include { KRAKENTOOLS_EXTRACT } from '../modules/ksumngs/nf-modules/krakentools/extract/main.nf'
include { KRAKENTOOLS_KREPORT2KRONA } from '../modules/ksumngs/nf-modules/krakentools/kreport2krona/main.nf'
include { KRONA_IMPORTTEXT } from '../modules/ksumngs/nf-modules/krona/importtext/main.nf'

workflow FILTERING {
    take:
    reads
    kraken2_db
    filter

    main:
    versions = Channel.empty()

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
        versions = versions.mix(KRAKENTOOLS_EXTRACT.out.versions)
    }

    KRAKENTOOLS_KREPORT2KRONA(KRAKEN2.out.kreport)
    KRONA_IMPORTTEXT(
        KRAKENTOOLS_KREPORT2KRONA.out.krona
            .map{ it.drop(1) }
            .collect()
    )
    KRONA_IMPORTTEXT.out.html.set{ krona }

    versions = versions.mix(KRAKEN2.out.versions)
    versions = versions.mix(KRONA_IMPORTTEXT.out.versions)
    versions = versions.mix(KRAKENTOOLS_KREPORT2KRONA.out.versions)

    emit:
    filtered
    log_out
    krona
    versions
}
