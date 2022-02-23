process RAXMLNG_PARSE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::raxml-ng=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:1.1.0--h32fcf60_0':
        'quay.io/biocontainers/raxml-ng:1.1.0--h32fcf60_0' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.raxml.rba"), emit: rba
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    raxml-ng \\
        --parse \\
        --threads auto{${task.cpus}} \\
        --msa ${alignment} \\
        --prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxmlng: \$(echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//')
    END_VERSIONS
    """
}
