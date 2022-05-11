#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CLIQUESNV_ILLUMINA } from '../../modules/ksumngs/nf-modules/cliquesnv/illumina/main.nf'
include { HAPLINK_HAPLOTYPES } from '../../modules/local/modules/haplink/haplotypes/main.nf'
include { HAPLINK_SEQUENCES } from '../../modules/local/modules/haplink/sequences/main.nf'
include { JSON2YAML } from '../../modules/local/json2yaml.nf'

workflow HAPLOTYPING {
    take:
    bam
    vcf
    reference

    main:
    versions = Channel.empty()

    if (params.platform == 'illumina') {
        CLIQUESNV_ILLUMINA(bam)
        CLIQUESNV_ILLUMINA.out.fasta.set{ fasta }
        versions = versions.mix(CLIQUESNV_ILLUMINA.out.versions)

        JSON2YAML(CLIQUESNV_ILLUMINA.out.json)
        JSON2YAML.out.yaml.set{ yaml }
        versions = versions.mix(JSON2YAML.out.versions)
    }
    else {
        HAPLINK_HAPLOTYPES(bam.join(vcf).join(reference))
        HAPLINK_HAPLOTYPES.out.yaml.set{ yaml }
        versions = versions.mix(HAPLINK_HAPLOTYPES.out.versions)

        HAPLINK_SEQUENCES(HAPLINK_HAPLOTYPES.out.yaml.join(reference))
        HAPLINK_SEQUENCES.out.fasta.set{ fasta }
        versions = versions.mix(HAPLINK_SEQUENCES.out.versions)
    }

    emit:
    yaml
    fasta
    versions
}
