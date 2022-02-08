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

/// summary: |
///   Collect some base stats on the quality of reads
/// input:
///   - tuple:
///       - name: sampleName
///         type: val(String)
///         description: Sample identifier
///       - name: reads
///         type: file
///         description: |
///           A single file (SE or interleaved) containing Illumina reads in fastq
///           format
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: FastQC report file
process FASTQC {
    label 'fastqc'

    input:
    tuple val(sampleName), file(reads)

    output:
    path("${sampleName}_fastqc.zip")

    script:
    """
    fastqc -t ${task.cpus} "${reads}"
    """
}

/// summary: |
///   Collect some base stats on the quality of reads
/// input:
///   - tuple:
///       - name: sampleName
///         type: val(String)
///         description: Sample identifier
///       - name: reads
///         type: file
///         description: File containing Nanopore reads in fastq format
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: NanoStat report file
process NANOSTAT {
    label 'nanostat'

    input:
    tuple val(sampleName), file(reads)

    output:
    path("${sampleName}_nanostat.log")

    script:
    """
    NanoStat -t ${task.cpus} --fastq "${reads}" > "${sampleName}_nanostat.log"
    """
}
