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
    tag "$accessionNumber"
    label 'run_local'
    label 'process_low'
    label 'error_backoff'

    conda (params.enable_conda ? 'bioconda::entrez-direct=16.2' : null)
    container 'docker.io/ncbi/edirect:12.5'

    // Prevent NCBI database timeouts by preventing this process from being run in
    // parallel. Trust me: this is actually faster than erroring out
    maxForks 1

    input:
    tuple val(meta), val(accessionNumber)

    output:
    tuple val(meta), path("${accessionNumber}.fasta"), emit: fasta
    tuple val(meta), path("${accessionNumber}.gb"), emit: genbank
    path "versions.yml", emit: versions

    script:
    """
    esearch \\
            -db nucleotide \\
            -query "${accessionNumber}" \\
        | efetch \\
            -format fasta \\
        > ${accessionNumber}.fasta
    esearch \\
            -db nucleotide \\
            -query "${accessionNumber}" \\
        | efetch \\
            -format gb \\
        > "${accessionNumber}.gb"

    # Fail on empty files, as efetch does not return failed exit codes for errors
    grep -q '[^[:space:]]' ${accessionNumber}.fasta || exit 1
    grep -q '[^[:space:]]' ${accessionNumber}.gb || exit 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
        efetch: \$(efetch -version)
    END_VERSIONS
    """
}
