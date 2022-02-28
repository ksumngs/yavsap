#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CAT_CAT } from '../modules/nf-core/modules/cat/cat/main.nf'
include { CAT_FASTQ } from '../modules/nf-core/modules/cat/fastq/main.nf'
include { MAFFT } from '../modules/nf-core/modules/mafft/main.nf'
include { RAXMLNG_BOOTSTRAP } from '../modules/ksumngs/nf-modules/raxmlng/bootstrap/main.nf'
include { RAXMLNG_PARSE } from '../modules/ksumngs/nf-modules/raxmlng/parse/main.nf'
include { RAXMLNG_SEARCH } from '../modules/ksumngs/nf-modules/raxmlng/search/main.nf'
include { RAXMLNG_SUPPORT } from '../modules/ksumngs/nf-modules/raxmlng/support/main.nf'
include { RENAME_HAPLOTYPES } from '../modules/local/rename-haplotypes.nf'
include { RENAME_NCBI } from '../modules/local/rename-ncbi.nf'

/// summary: Create a phylogenetic tree
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: Sample identifier
///       - name: Alignment
///         type: file
///         description: |
///           Plain-text multi-alignment file. RAxML-NG supports FASTA, PHYLIP, and
///           CATG formats
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: Annotated support tree
workflow PHYLOGENETIC_TREE {
    take:
    sequences
    consensus
    genomes
    genome_table

    main:
    versions = Channel.empty()

    RENAME_NCBI(genomes, genome_table)

    CAT_FASTQ(
        sequences
        .join(consensus)
        .map{ [
            ['id':it[0].id, 'single_end': true, 'strandedness': null],
            [it[1], it[2]]
        ] }
    )

    RENAME_HAPLOTYPES(CAT_FASTQ.out.reads)

    CAT_CAT(
        RENAME_NCBI.out.fasta
            .mix(RENAME_HAPLOTYPES.out.fasta.map{ it.drop(1) })
            .collect(),
        'genotypes.fasta'
    )

    MAFFT(
        CAT_CAT.out.file_out.map{ [
            ['id': 'collective', 'single_end': null, 'strandedness': null],
            it
        ] }
    )

    RAXMLNG_PARSE(MAFFT.out.fas.map{ it[1] })

    RAXMLNG_SEARCH(RAXMLNG_PARSE.out.rba)
    RAXMLNG_BOOTSTRAP(RAXMLNG_PARSE.out.rba)

    RAXMLNG_SUPPORT(RAXMLNG_SEARCH.out.best_tree, RAXMLNG_BOOTSTRAP.out.bootstraps)

    tree = RAXMLNG_SUPPORT.out.support

    versions = versions.mix(RENAME_NCBI.out.versions)
    versions = versions.mix(CAT_FASTQ.out.versions)
    versions = versions.mix(RENAME_HAPLOTYPES.out.versions)
    versions = versions.mix(CAT_CAT.out.versions)
    versions = versions.mix(MAFFT.out.versions)
    versions = versions.mix(RAXMLNG_PARSE.out.versions)
    versions = versions.mix(RAXMLNG_SEARCH.out.versions)
    versions = versions.mix(RAXMLNG_BOOTSTRAP.out.versions)
    versions = versions.mix(RAXMLNG_SUPPORT.out.versions)

    emit:
    tree
    versions
}
