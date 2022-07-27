process GFFCAT {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/millironx/biojulia:1.6.7-2.0.5-387d929'

    input:
    tuple val(meta), path(gffs, stageAs: "input**.gff")

    output:
    tuple val(meta), path("*.gff"), emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gffcat \\
        '${meta.accession_num}' \\
        ${gffs} \\
        > ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        julia: \$(julia -v | awk '{print \$3}')
        "GFF3.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("00701ae9-d1dc-5365-b64a-a3a3ebf5695e")].version))')
        "HTTP.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("af1dc308-cb6b-11e8-32f0-31192efa90f6")].version))')
    END_VERSIONS
    """
}
