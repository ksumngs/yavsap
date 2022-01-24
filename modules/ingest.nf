#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/// summary: |
///   Takes a collection of single-end sequencing reads and converts them into a
///   single file with a safe sample and filename
/// input:
///   - tuple:
///       - name: givenName
///         type: val(String)
///         description: Identifier for this sample
///       - name: reads
///         type: file
///         description: Collection of individual reads files
/// output:
///   - tuple:
///       - type: val(String)
///         description: Identifier for this sample, cleaned of any escape characters
///       - type: path
///         description: The single, complete reads file
process SINGLE_PREPROCESS {
    label 'parallelzip'

    input:
    tuple val(givenName), file(reads)

    output:
    tuple val(sampleName), path("${sampleName}.fastq.gz")

    script:
    sampleName = givenName.split('_')[0]
    """
    parallel -j${task.cpus} gunzip -f ::: *.gz
    cat *.f*q > "${sampleName}.fastq"
    pigz -p${task.cpus} "${sampleName}.fastq"
    """
}
