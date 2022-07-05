process EDIRECT_EFETCH {
    tag "$search"
    label 'process_low'
    label 'run_local'
    label 'error_backoff'

    conda (params.enable_conda ? "bioconda::entrez-direct=16.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1' :
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"

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
