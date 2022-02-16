process XTRACT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::entrez-direct=16.2' : null)
    container 'docker.io/ncbi/edirect:12.5'

    // Prevent NCBI database timeouts by preventing this process from being run in
    // parallel. Trust me: this is actually faster than erroring out
    maxForks 1

    input:
    tuple val(meta), path(xml), val(pattern), val(element)

    output:
    tuple val(meta), env(RESULT), emit: result
    path "versions.yml"        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export RESULT=$(\\
        xtract \\
            < ${xml} \\
            -pattern ${pattern} \\
            -element ${element} \\
            ${args} \\
        )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xtract: \$(xtract -version)
    END_VERSIONS
    """
}
