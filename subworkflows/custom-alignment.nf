include { MINIMAP2_ALIGN } from '../modules/ksumngs/nf-modules/minimap2/align/main.nf'
include { SAMTOOLS_SORT } from '../modules/nf-core/modules/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main.nf'

workflow CUSTOM_ALIGNMENT {
    take:
    reads

    main:
    versions = Channel.empty()

    // Realign reads to the reference genome
    // Note: Normally, minimap2 outputs paf, but we have forced it to output sam via
    // ext.args in modules.config
    MINIMAP2_ALIGN(reads)
    SAMTOOLS_SORT(MINIMAP2_ALIGN.out.paf)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    SAMTOOLS_SORT.out.bam.set{ bam }
    SAMTOOLS_INDEX.out.bai.set{ bai }

    versions = versions.mix(MINIMAP2_ALIGN.out.versions)
    versions = versions.mix(SAMTOOLS_SORT.out.versions)
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)


    emit:
    bam
    bai
    versions
}
