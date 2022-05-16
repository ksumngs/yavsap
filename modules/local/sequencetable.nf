process SEQUENCETABLE {
    tag "$haplotypes"
    label 'process_low'
    cache false

    container 'quay.io/millironx/biojulia:1.6.6-2.0.5-9877308'

    input:
    path(haplotypes, stageAs: 'haplotypes.yml')
    path(reference, stageAs: 'reference.fasta')
    path(template, stageAs: 'template.html')
    path(freezetable_js, stageAs: 'freezetable.jquery.js')

    output:
    path "*_mqc.html", emit: mqc_html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sequencetable \\
            ${haplotypes} \\
            ${reference} \\
            ${template} \\
            ${freezetable_js} \\
            sequencetable_mqc.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        julia: \$(julia -v | awk '{print \$3}')
        "EzXML.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615")].version))')
        "FASTX.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("c2308a5c-f048-11e8-3e8a-31650f418d12")].version))')
        "Kelpie.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("1b112299-d6bc-44e2-912a-478f25731460")].version))')
        "Mustache.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("ffc61752-8dc7-55ee-8c37-f3e9cdd09e70")].version))')
        "YAML.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("ddb6d928-2868-570f-bddf-ab3f9cf99eb6")].version))')
    END_VERSIONS
    """
}
