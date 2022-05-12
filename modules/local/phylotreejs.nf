process PHYLOTREEJS {
    tag "${newick}"
    label 'process_low'

    container 'quay.io/millironx/juliapro:1.6.6-2ed3693'

    input:
    path(newick, stageAs: 'tree.nwk')
    path(template, stageAs: 'template.html')
    path(css, stageAs: 'styles.css')
    path(d3, stageAs: 'd3.js')
    path(underscore, stageAs: 'underscore.js')
    path(phylotree, stageAs: 'phylotree.js')

    output:
    path "*_mqc.html", emit: mqc_html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    phylotreegen \\
        ${newick} \\
        ${template} \\
        phylotree_mqc.html \\
        ${css} \\
        ${d3} \\
        ${underscore} \\
        ${phylotree}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        julia: \$(julia -v | awk '{print \$3}')
    END_VERSIONS
    """
}
