process HAPLINK_SEQUENCES {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/millironx/haplink:0.5.1'

    input:
    tuple val(meta), file(yaml), file(reference)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplink sequences \\
        --haplotypes ${yaml} \\
        --reference ${reference} \\
        --output ${prefix}.haplotypes.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplink: \$(haplink --version)
    END_VERSIONS
    """
}
