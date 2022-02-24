process RAXMLNG_SUPPORT {
    tag "$tree"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::raxml-ng=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:1.1.0--h32fcf60_0':
        'quay.io/biocontainers/raxml-ng:1.1.0--h32fcf60_0' }"

    input:
    path tree
    path bootstraps

    output:
    path "*.raxml.support", emit: bam
    path "versions.yml"   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    raxml-ng \\
        --support \\
        --threads auto{${task.cpus}} \\
        --tree ${tree} \\
        --bs-trees ${bootstraps} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxmlng: \$(echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//')
    END_VERSIONS
    """
}
