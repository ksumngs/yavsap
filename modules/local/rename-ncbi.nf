process RENAME_NCBI {
    tag "$fasta"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::r-phylotools=0.2.2" : null)
    container 'docker.io/biocontainers/phylotools:v0.2.4_cv1'

    input:
    path(fasta)
    path(tsv)

    output:
    path "renamed.fasta", emit: fasta
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    renamerefseqs ${fasta} ${tsv} renamed.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
