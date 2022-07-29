process HAPLINK_CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/millironx/haplink:0.7.1'

    input:
    tuple val(meta), path(variantcalls), path(reference)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplink consensus \\
        --reference ${reference} \\
        --variants ${variantcalls} \\
        --prefix ${prefix} \\
        --output ${prefix}.consensus.fasta \\
        ${args} --julia-args -t${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplink: \$(haplink --version)
    END_VERSIONS
    """
}
