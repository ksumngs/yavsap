#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/// summary: |
///   Downloads a genome in fasta and genbank format from the NCBI databases
/// input:
///   - name: accessionNumber
///     type: val(String)
///     description: The accession number or gid to download from NCBI
/// output:
///   - name: fasta
///     type: path
///     description: The requested genome in fasta format
///   - name: genbank
///     type: path
///     description: The requested genome in GenBank format
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
