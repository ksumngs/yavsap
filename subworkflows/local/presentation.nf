include { IGV } from '../../modules/local/igv'
include { PHYLOTREEJS } from '../../modules/local/phylotreejs'

workflow PRESENTATION {
    take:
    bam
    tree

    main:
    versions = Channel.empty()

    igv_js = file(params.igv_js, checkIfExists: true)
    igv_template = file("${workflow.projectDir}/assets/igv_mqc.html", checkIfExists: true)
    IGV(
        bam
            .map{ "${it[0].id}" }
            .collectFile(name: 'samplenames.txt', newLine: true),
        igv_js,
        igv_template
    )
    IGV.out.mqc_html.set{ igv }
    versions = versions.mix(IGV.out.versions)

    phylotree_css = file(params.phylotree_css, checkIfExists: true)
    d3_js = file(params.d3_js, checkIfExists: true)
    underscore_js = file(params.underscore_js, checkIfExists: true)
    phylotree_js = file(params.phylotree_js, checkIfExists: true)
    phylotree_template = file(
        "${workflow.projectDir}/assets/phylotree_mqc.html", checkIfExists: true
    )
    PHYLOTREEJS(
        tree, phylotree_template, phylotree_css, d3_js, underscore_js, phylotree_js
    )
    PHYLOTREEJS.out.mqc_html.set{ phylotree }
    versions = versions.mix(PHYLOTREEJS.out.versions)

    emit:
    igv
    phylotree
    versions
}
