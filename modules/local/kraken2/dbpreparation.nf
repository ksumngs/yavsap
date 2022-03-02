process KRAKEN2_DBPREPARATION {
    tag "$db"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pyaml=15.8.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyaml:15.8.2--py36_0' :
        'quay.io/biocontainers/pyaml:15.8.2--py36_0' }"

    input:
    path db

    output:
    path "${dbname}", emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    dbname = db.getSimpleName()
    """
    find . -name "*.tar.gz" -exec tar -xvf {} \\;
    find . -name "*.tgz" -exec tar -xvf {} \\;
    mkdir ${dbname}
    find . -name "*.k2d" -exec mv {} ${dbname} \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
