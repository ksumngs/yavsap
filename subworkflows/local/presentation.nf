include { HAPLOTYPECONVERT } from '../../modules/local/haplotypeconvert'
include { IGV } from '../../modules/local/igv'
include { PHYLOTREEJS } from '../../modules/local/phylotreejs'
include { SEQUENCETABLE } from '../../modules/local/sequencetable'

workflow PRESENTATION {
    take:
    bam
    reference
    strain
    accession
    consensus
    haplotype_fasta
    haplotype_yaml
    tree

    main:
    versions = Channel.empty()

    HAPLOTYPECONVERT(
        strain
            .join(accession)
            .join(consensus)
            .join(haplotype_fasta, remainder: true)
            .join(haplotype_yaml, remainder: true),
        reference
    )
    HAPLOTYPECONVERT
        .out
        .yaml
        .map{ it[1] }
        .collectFile(name: 'collated_haplotypes.yml', newLine: true)
        .set{ ch_collected_haplotypes }
    versions = versions.mix(HAPLOTYPECONVERT.out.versions)

    freezetable_js = file(params.freezetable_js, checkIfExists: true)
    sequencetable_template = file(
        "${workflow.projectDir}/assets/kelpie_mqc.html", checkIfExists: true
    )
    SEQUENCETABLE(
        ch_collected_haplotypes,
        reference,
        sequencetable_template,
        freezetable_js
    )
    SEQUENCETABLE.out.mqc_html.set{ seqtable }
    versions = versions.mix(SEQUENCETABLE.out.versions)

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
    seqtable
    igv
    phylotree
    versions
}
