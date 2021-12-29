#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/// summary: |
///   Downloads
process NCBI_DOWNLOAD {
    label 'edirect'
    label 'run_local'
    label 'process_low'
    label 'error_backoff'

    input:
    val(accessionNumber)

    output:
    path("${accessionNumber}.fasta"), emit: fasta
    path("${accessionNumber}.gb"), emit: genbank

    script:
    """
    efetch \
            -db nucleotide \
            -id ${accessionNumber} \
            -format fasta \
        > ${accessionNumber}.fasta

    # Avoid NCBI timeouts
    sleep 0.3

    efetch \
            -db nucleotide \
            -id ${accessionNumber} \
            -format gb \
        > ${accessionNumber}.gb

    # Fail on empty files, as efetch does not return failed exit codes for errors
    grep -q '[^[:space:]]' ${accessionNumber}.fasta || exit 1
    grep -q '[^[:space:]]' ${accessionNumber}.gb || exit 1
    """
}
