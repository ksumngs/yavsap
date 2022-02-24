process CLIQUESNV_ILLUMINA {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::cliquesnv=2.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cliquesnv:2.0.3--hdfd78af_0':
        'quay.io/biocontainers/cliquesnv:2.0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.json") , emit: json
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def jmemstring = task.memory.toMega() + 'M'
    """
    cliquesnv \\
        -Xmx${jmemstring} \\
        -threads ${task.cpus} \\
        -m snv-illumina \\
        -in ${bam} \\
        ${args} \\
        -outDir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cliquesnv: \$(cliquesnv | head -n1 | sed 's/CliqueSNV version: //')
    END_VERSIONS
    """
}
