process RENAME_HAPLOTYPES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::r-phylotools=0.2.2" : null)
    container 'docker.io/biocontainers/phylotools:v0.2.4_cv1'

    input:
    tuple val(meta), path(fasta, stageAs: 'input.fasta')

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    renamehapseqs input.fasta ${prefix} ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
