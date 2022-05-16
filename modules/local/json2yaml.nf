process JSON2YAML {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pyaml=15.8.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyaml:15.8.2--py36_0' :
        'quay.io/biocontainers/pyaml:15.8.2--py36_0' }"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*.yaml"), emit: yaml
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    json2yaml ${json} > ${prefix}.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g' | sed 's/ ::.*//g')
    END_VERSIONS
    """
}
