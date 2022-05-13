process HAPLOTYPECONVERT {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/millironx/biojulia:1.6.6-2.0.5-9877308'

    input:
    tuple val(meta), val(strain), val(ncbi), path(consensus), path(haplotype_fasta), path(haplotype_yaml)
    path(reference)

    output:
    tuple val(meta), path("*.yaml"), emit: yaml
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta = haplotype_fasta ?: 'nothing.fasta'
    def yaml = haplotype_yaml ?: 'nothing.yaml'
    """
    haplotypestandardizer \\
        ${meta.id} \\
        ${reference} \\
        ${strain} \\
        ${ncbi} \\
        ${consensus} \\
        ${yaml} \\
        ${fasta} \\
        ${prefix}.haplotypes.std.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        julia: \$(julia -v | awk '{print \$3}')
    END_VERSIONS
    """
}
