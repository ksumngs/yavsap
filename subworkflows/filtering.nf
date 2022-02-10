#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

import static java.lang.Math.sqrt
import static java.lang.Math.round

KrakenDbSize = file(params.kraken2_db).toFile().directorySize()
KrakenAllocationSize = round(sqrt(KrakenDbSize) + KrakenDbSize)

workflow read_filtering {
    take:
    InputReads

    main:
    if (params.skip_filtering) {
        FilteredReads = InputReads
        KrakenReports = Channel.from([])
    }
    else {
        classification(InputReads)
        KrakenReports = classification.out

        KrakenReads = InputReads.join(KrakenReports)

        // Filter out the non-viral reads
        filtering(KrakenReads)
        FilteredReads = filtering.out
    }

    emit:
    FilteredReads
    KrakenReports
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
