include { CLIQUESNV_ILLUMINAVC } from '../../modules/ksumngs/nf-modules/cliquesnv/illuminavc/main.nf'
include { HAPLINK_VARIANTS } from '../../modules/local/modules/haplink/variants/main.nf'

workflow VARIANTS {
    take:
    bam
    bai
    reference

    main:
    versions = Channel.empty()
    vcf = Channel.empty()

    if (params.platform == 'illumina') {
        CLIQUESNV_ILLUMINAVC(bam)
        CLIQUESNV_ILLUMINAVC.out.vcf.set{ vcf }
        versions = versions.mix(CLIQUESNV_ILLUMINAVC.out.versions.first())
    }

    if (params.platform == 'nanopore') {
        HAPLINK_VARIANTS(bam.join(bai).join(reference))
        HAPLINK_VARIANTS.out.vcf.set{ vcf }
        versions = versions.mix(HAPLINK_VARIANTS.out.versions.first())
    }

    emit:
    vcf
    versions
}
