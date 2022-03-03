process RAXMLNG_SEARCH {
    tag "$msa"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::raxml-ng=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:1.1.0--h32fcf60_0':
        'quay.io/biocontainers/raxml-ng:1.1.0--h32fcf60_0' }"

    input:
    path msa

    output:
    path "*.raxml.bestTree", emit: best_tree
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    raxml-ng \\
        --threads auto{${task.cpus}} \\
        --workers auto \\
        --msa ${msa} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxmlng: \$(echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//')
    END_VERSIONS
    """
}
