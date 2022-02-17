process NANOSTAT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::nanostat=1.6.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/nanostat:1.6.0--pyhdfd78af_0' :
    'quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    NanoStat \\
            -t ${task.cpus} \\
            --fastq ${reads} \\
            ${args} \\
        > ${prefix}_nanostat.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanostat: \$(NanoStat -v | sed 's/NanoStat //')
    END_VERSIONS
    """
}
