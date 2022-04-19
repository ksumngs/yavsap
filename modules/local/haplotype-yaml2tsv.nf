process HAPLOTYPE_YAML2TSV {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pyaml=15.8.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyaml:15.8.2--py36_0' :
        'quay.io/biocontainers/pyaml:15.8.2--py36_0' }"

    input:
    tuple val(meta), path(yaml)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplotype-parser ${prefix} ${yaml} ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
