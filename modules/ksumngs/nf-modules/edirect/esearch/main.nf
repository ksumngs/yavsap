process EDIRECT_ESEARCH {
    tag "${ (query.length() > 100) ? query.substring(0, 100) + '...' : query }"
    label 'process_low'
    label 'run_local'
    label 'error_backoff'

    conda (params.enable_conda ? "bioconda::entrez-direct=16.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1' :
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"

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
