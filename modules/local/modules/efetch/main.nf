process EFETCH {
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
    tuple val(meta), path(search), val(format), val(mode)

    output:
    tuple val(meta), path("*.${mode ?: format ?: 'txt'}"), emit: txt
    path "versions.yml"                                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = mode ?: format ?: 'txt'
    def format_flag = format ? "-format ${format}" : ''
    def mode_flag = mode ? "-mode ${mode}" : ''
    """
    efetch \\
            < ${search} \\
            ${format_flag} \\
            ${mode_flag} \\
            ${args} \\
        > ${prefix}.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        efetch: \$(efetch -version)
    END_VERSIONS
    """
}
