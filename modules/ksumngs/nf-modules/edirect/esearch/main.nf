process EDIRECT_ESEARCH {
    tag "$query"
    label 'run_local'
    label 'process_low'
    label 'error_backoff'

    conda (params.enable_conda ? "bioconda::entrez-direct=16.2" : null)
    container 'docker.io/ncbi/edirect:12.5'

    // Prevent NCBI database timeouts by preventing this process from being run in
    // parallel. Trust me: this is actually faster than erroring out
    maxForks 1

    input:
    val(query)
    val(db)

    output:
    path "search.xml"   , emit: xml
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    esearch \\
            -db ${db} \\
            -query "${query}" \\
            ${args} \\
        > search.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
    END_VERSIONS
    """
}
