process HAPLINK_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/millironx/haplink:0.7.1'

    input:
    tuple val(meta), file(bam), file(bai), file(reference)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplink variants \\
        --bam ${bam} \\
        --reference ${reference} \\
        --output ${prefix}.vcf \\
        ${args} --julia-args -t${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplink: \$(haplink --version)
    END_VERSIONS
    """
}
