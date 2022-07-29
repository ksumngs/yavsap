process HAPLINK_HAPLOTYPES {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/millironx/haplink:0.7.1'

    input:
    tuple val(meta), path(bam), path(vcf), path(reference)

    output:
    tuple val(meta), path("*.yaml"), emit: yaml
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplink haplotypes \\
        --bam "${bam}" \\
        --variants "${vcf}" \\
        --reference "${reference}" \\
        --output "${prefix}.haplotypes.yaml" \\
        ${args} --julia-args -t${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplink: \$(haplink --version)
    END_VERSIONS
    """
}
