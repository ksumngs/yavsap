process ESEARCH {
    tag "$meta.id"
    label 'run_local'
    label 'process_low'
    label 'error_backoff'

    conda (params.enable_conda ? 'bioconda::entrez-direct=16.2' : null)
    container 'docker.io/ncbi/edirect:12.5'

    // Prevent NCBI database timeouts by preventing this process from being run in
    // parallel. Trust me: this is actually faster than erroring out
    maxForks 1

    input:
    tuple val(meta), val(db), val(query)

    output:
    tuple val(meta), path("*.xml"), emit: xml
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    esearch \\
            -db "${db}" \\
            -query "${query}" \\
            ${args} \\
        > "${prefix}.xml"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
    END_VERSIONS
    """
}
