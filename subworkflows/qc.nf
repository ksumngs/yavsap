#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/// summary: |
///   Merge paired-end fastq reads files into a single file
/// input:
///   - tuple:
///       - name: sampleName
///         type: val(String)
///         description: Sample identifier
///       - name: reads
///         type: file
///         description: Set of paired-end reads files
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: Interleaved paired-end reads file
process SEQTK_INTERLEAVE {
    label 'seqtk'
    label 'process_low'

    input:
    tuple val(sampleName), file(reads)

    output:
    tuple val(sampleName), path("${sampleName}.fastq.gz")

    script:
    """
    seqtk mergepe "${reads}" | gzip > "${sampleName}.fastq.gz"
    """
}
