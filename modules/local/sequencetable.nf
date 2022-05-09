process SEQUENCETABLE {
    tag "$tsv"
    label 'process_low'
    cache false

    container 'quay.io/millironx/biojulia:1.6.6-1.1.4-9409225'

    input:
    path tsv
    path sam
    path reference
    path tree
    path multiqc
    path krona

    output:
    path "*.html", emit: html
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    multiqc_flag = multiqc ? '--multiqc' : '--no-multiqc'
    krona_flag = krona ? '--krona' : '--no-krona'
    tree_flag = tree ? "--newick ${tree}" : ''
    """
    sequence-table \\
            ${tsv} \\
            ${sam} \\
            ${reference} \\
            ${multiqc_flag} \\
            ${krona_flag} \\
            ${tree_flag} \\
        > seq-table.partial.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        julia: \$(julia -v | awk '{print \$3}')
    END_VERSIONS
    """
}
