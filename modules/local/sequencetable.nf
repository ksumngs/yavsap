process SEQUENCETABLE {
    tag "$tsv"
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
    END_VERSIONS
    """
}
