process IGV {
    tag "$samplelist"
    label 'process_low'

    container 'quay.io/millironx/juliapro:1.6.6-2ed3693'
    input:
    path(samplelist, stageAs: 'samples.txt')
    path(igvjs, stageAs: 'igv.js')
    path(template, stageAs: 'template.html')

    output:
    path "*_mqc.html", emit: mqc_html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    igvgen \\
        ${samplelist} \\
        ${igvjs} \\
        ${template} \\
        igv_mqc.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        julia: \$(julia -v | awk '{print \$3}')
    END_VERSIONS
    """
}
