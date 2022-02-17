include { MINIMAP2_ALIGN } from '../modules/nf-core/modules/minimap2/align/main.nf'
include { SAMTOOLS_SORT } from '../modules/nf-core/modules/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main.nf'

workflow ALIGNMENT {
    take:
    reads
    reference

    main:
    // Realign reads to the reference genome
    // Note: Normally, minimap2 outputs paf, but we have forced it to output sam via
    // ext.args in modules.config
    MINIMAP2_ALIGN(reads, reference)
    SAMTOOLS_SORT(MINIMAP2_ALIGN.out.paf)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    SAMTOOLS_SORT.out.bam.set{ bam }
    SAMTOOLS_INDEX.out.bai.set{ bai }

    emit:
    bam
    bai
}
