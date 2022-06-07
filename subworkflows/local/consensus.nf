include { CLIQUESNV_CONSENSUSILLUMINA } from '../../modules/ksumngs/nf-modules/cliquesnv/consensusillumina/main.nf'
include { HAPLINK_CONSENSUS } from '../../modules/local/haplink/consensus'
include { HAPLINK_VARIANTS } from '../../modules/local/haplink/variants'

workflow CONSENSUS {
    take:
    bam
    bai
    reference

    main:
    versions = Channel.empty()

    if (params.platform == 'illumina') {
        CLIQUESNV_CONSENSUSILLUMINA(bam)
        CLIQUESNV_CONSENSUSILLUMINA.out.fasta.set{ fasta }
        versions = versions.mix(CLIQUESNV_CONSENSUSILLUMINA.out.versions)
    }
    else if (params.platform == 'nanopore') {
        BamPlusReference = bam
            .join(bai)
            .combine(reference)

        HAPLINK_VARIANTS(BamPlusReference)
        HAPLINK_VARIANTS.out.vcf.set{ vcf }
        versions = versions.mix(HAPLINK_VARIANTS.out.versions)

        HAPLINK_CONSENSUS(vcf.combine(reference))
        HAPLINK_CONSENSUS.out.fasta.set{ fasta }
        versions = versions.mix(HAPLINK_CONSENSUS.out.versions)
    }

    emit:
    fasta
    versions
}
