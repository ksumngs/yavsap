#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CLIQUESNV_ILLUMINA } from '../modules/ksumngs/nf-modules/cliquesnv/illumina/main.nf'
include { CLIQUESNV_ILLUMINAVC } from '../modules/ksumngs/nf-modules/cliquesnv/illuminavc/main.nf'
include { HAPLINK_HAPLOTYPES } from '../modules/local/modules/haplink/haplotypes/main.nf'
include { HAPLINK_SEQUENCES } from '../modules/local/modules/haplink/sequences/main.nf'
include { HAPLINK_VARIANTS } from '../modules/local/modules/haplink/variants/main.nf'
include { JSON2YAML } from '../modules/local/json2yaml.nf'
include { PHYLOGENETIC_TREE } from './phylogenetics.nf'

workflow HAPLOTYPING {
    take:
    alignments
    references

    main:
    if (params.platform == 'illumina') {
        // Drop the BAM index: CliqueSNV doesn't need it
        alignments
            .map{ it.dropRight(1) }
            .set{ UnindexedAlignments }

        // Do variant calling
        CLIQUESNV_ILLUMINAVC(UnindexedAlignments)
        CLIQUESNV_ILLUMINAVC.out.vcf.set{ vcf }

        // Do haplotype calling
        CLIQUESNV_ILLUMINA(UnindexedAlignments)
        CLIQUESNV_ILLUMINA.out.fasta.set{ fasta }

        // Convert haplotyp JSON to YAML
        JSON2YAML(CLIQUESNV_ILLUMINA.out.json)
        JSON2YAML.out.yaml.set{ yaml }
    }
    else {
        HAPLINK_VARIANTS(alignments.join(references))
        HAPLINK_VARIANTS.out.vcf.set{ vcf }

        HAPLINK_HAPLOTYPES(
            alignments
                .map{ it.dropRight(1) }
                .join(HAPLINK_VARIANTS.out.vcf)
        )
        HAPLINK_HAPLOTYPES.out.yaml.set{ yaml }

        HAPLINK_SEQUENCES(
            HAPLINK_HAPLOTYPES.out.yaml
                .join(references)
        )
        HAPLINK_SEQUENCES.out.fasta.set{ fasta }
    }

    // AllHapSequences = HaplotypeSequences.join(ConsensusSequences)

    // alignment(AllHapSequences, StrainGenomes) | \
    //     PHYLOGENETIC_TREE

    // trees = PHYLOGENETIC_TREE.out

    emit:
    vcf
    yaml
    fasta
}
