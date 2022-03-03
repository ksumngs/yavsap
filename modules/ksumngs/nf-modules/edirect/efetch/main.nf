process EDIRECT_EFETCH {
    tag "$search"
    label 'run_local'
    label 'process_low'
    label 'error_backoff'

    conda (params.enable_conda ? "bioconda::entrez-direct=16.2" : null)
    container 'docker.io/ncbi/edirect:12.5'

    input:
    path(search)
    val(format)
    val(mode)

    output:
    path "*.${mode ?: format ?: 'txt'}", emit: txt
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ext = mode ?: format ?: 'txt'
    def format_flag = format ? "-format ${format}" : ''
    def mode_flag = mode ? "-mode ${mode}" : ''
    """
    efetch \\
            < ${search} \\
            ${format_flag} \\
            ${mode_flag} \\
            ${args} \\
        > result.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        efetch: \$(efetch -version)
    END_VERSIONS
    """
}
