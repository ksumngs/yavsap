process VARTABLE {
    tag "${metadata}"
    label 'process_low'

    container 'quay.io/millironx/biojulia:1.6.7-2.0.5-bb3c4be'

    input:
    path(variants)
    path(metadata, stageAs: 'tool.yml')
    path(template, stageAs: 'template.html')

    output:
    path "*_mqc.html", emit: mqc_html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    vartable \\
        ${metadata} \\
        ${template} \\
        ${variants} \\
        > vartable_mqc.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        julia: \$(julia -v | awk '{print \$3}')
        "EzXML.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615")].version))')
        "Kelpie.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("1b112299-d6bc-44e2-912a-478f25731460")].version))')
        "Mustache.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("ffc61752-8dc7-55ee-8c37-f3e9cdd09e70")].version))')
        "OrderedCollections.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("bac558e1-5e72-5ebc-8fee-abe8a469f55d")].version))')
        "VariantCallFormat.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("28eba6e3-a997-4ad9-87c6-d933b8bca6c1")].version))')
        "YAML.jl": \$(julia -e 'using Pkg, UUIDs; println(string(Pkg.dependencies()[UUID("ddb6d928-2868-570f-bddf-ab3f9cf99eb6")].version))')
    END_VERSIONS
    """
}
