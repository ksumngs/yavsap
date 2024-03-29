include { HAPLOTYPECONVERT } from '../../modules/local/haplotypeconvert'
include { IGV } from '../../modules/local/igv'
include { PHYLOTREEJS } from '../../modules/local/phylotreejs'
include { SEQUENCETABLE } from '../../modules/local/sequencetable'
include { VARTABLE } from '../../modules/local/vartable'

workflow PRESENTATION {
    take:
    bam
    reference
    gff
    strain
    consensus
    variants
    haplotype_fasta
    haplotype_yaml
    tree

    main:
    versions = Channel.empty()

    tool_meta = []
    if (params.platform == 'illumina') {
        tool_meta = file(
            "${workflow.projectDir}/assets/cliquesnv_info.yml", checkIfExists: true
        )
    }
    else if (params.platform == 'nanopore') {
        tool_meta = file(
            "${workflow.projectDir}/assets/haplink_info.yml", checkIfExists: true
        )
    }
    vartable_template = file(
        "${workflow.projectDir}/assets/variants_mqc.html", checkIfExists: true
    )

    VARTABLE(variants.collect{ it[1] }, tool_meta, vartable_template)
    VARTABLE.out.mqc_html.set{ vartable }
    versions = versions.mix(VARTABLE.out.versions)

    strain
        .map{ [it[0].id, it] }
        .join(
            consensus.map{ [it[0].id, it] }
        )
        .map{ [it[1][0], it[2][1]] }
        .join(haplotype_fasta)
        .join(haplotype_yaml)
        .map{ it[2] ? it : [it[0], it[1], [], it[3]] }
        .map{ it[3] ? it : [it[0], it[1], it[2], it[3]] }
        .set{ ch_haplotype_meta }
    ch_haplotype_meta.dump(tag: 'haplotype_meta')

    HAPLOTYPECONVERT(ch_haplotype_meta, reference.map{ it[1] })
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
        reference.map{ it[1] },
        gff.map{ it[1] },
        sequencetable_template,
        tool_meta,
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
    vartable
    seqtable
    igv
    phylotree
    versions
}
