#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

import static java.lang.Math.sqrt
import static java.lang.Math.round

KrakenDbSize = file(params.kraken2_db).toFile().directorySize()
KrakenAllocationSize = round(sqrt(KrakenDbSize) + KrakenDbSize)

include { KRAKEN2 } from '../modules/local/modules/kraken2/main.nf'
include { KRAKENTOOLS_EXTRACT } from '../modules/local/modules/krakentools/extract/main.nf'

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

// Classify reads using Kraken
process classification {
    label 'kraken'
    label 'process_high_memory'
    memory "${KrakenAllocationSize} B"
    publishDir "${params.outdir}/classification", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(readsFile)

    output:
    tuple val(sampleName), file("${sampleName}.kraken"), file("${sampleName}.kreport")

    script:
    pairedflag = params.paired ? '--paired' : ''
    """
    kraken2 --db ${params.kraken2_db} --threads ${task.cpus} \
        --report "${sampleName}.kreport" \
        --output "${sampleName}.kraken" \
        ${pairedflag} \
        ${readsFile}
    """
}

// Pull the viral reads and any unclassified reads from the original reads
// files for futher downstream processing using KrakenTools
process filtering {
    label 'krakentools'
    label 'process_low'

    input:
    tuple val(sampleName), file(readsFile), file(krakenFile), file(krakenReport)

    output:
    tuple val(sampleName), file("${sampleName}_filtered*.fastq.gz")

    script:
    read2flagin  = (params.paired) ? "-s2 ${readsFile[1]}" : ''
    read1flagout = (params.paired) ? "-o ${sampleName}_filtered_R1.fastq" : " -o ${sampleName}_filtered.fastq"
    read2flagout = (params.paired) ? "-o2 ${sampleName}_filtered_R2.fastq" : ''
    """
    extract_kraken_reads.py -k ${krakenFile} \
        -s ${readsFile[0]} ${read2flagin} \
        -r ${krakenReport} \
        -t ${params.keep_taxid} --include-children \
        --fastq-output \
        ${read1flagout} ${read2flagout}
    gzip ${sampleName}_filtered*.fastq
    """
}
