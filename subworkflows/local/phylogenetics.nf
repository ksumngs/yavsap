#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MAFFT } from '../../modules/nf-core/modules/mafft/main.nf'
include { RAXMLNG_BOOTSTRAP } from '../../modules/ksumngs/nf-modules/raxmlng/bootstrap/main.nf'
include { RAXMLNG_PARSE } from '../../modules/ksumngs/nf-modules/raxmlng/parse/main.nf'
include { RAXMLNG_SEARCH } from '../../modules/ksumngs/nf-modules/raxmlng/search/main.nf'
include { RAXMLNG_SUPPORT } from '../../modules/ksumngs/nf-modules/raxmlng/support/main.nf'

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
    haplotype_fasta
    consensus_fasta
    genome_fasta
    genome_strain

    main:
    versions = Channel.empty()

    // This beautiful one-liner deserves docs of a full-blown process.
    // In fact, it actually cuts out the need for a class of finicky R
    // processes, so here it is step-by-step
    // 1. Split the id'ed fasta haplotypes into record objects
    // 2. Remove any consensus sequences (more of a concern with HapLink)
    // 3. Reassign an id to HapLink haplotypes based on the SHA1 hash that
    //    HapLink assigns to each haplotype
    // 4. Reasign an id to CliqueSNV haplotypes based on the ID_NUM_FREQ
    //    extended tag that CliqueSNV assigns to each haplotype
    // 5. Create a fasta string with all the sequences
    // 6. Create a file out of it
    haplotype_fasta
        .splitFasta(record: [id: true, sequence: true], elem: 1)
        .filter{ !it[1].id.toLowerCase().contains('consensus') }
        .map{ [it[0], [id:it[1].id, sequence:it[1].sequence]] }
        .map{
            [
                it[0],
                [
                    id: it[1].id ==~ /^[0-9a-f]{8}$/ ?
                        "${it[0].id}_haplotype_${it[1].id}" : it[1].id,
                    sequence: it[1].sequence
                ]
            ]
        }
        .map{
            [
                it[0],
                [
                    id: it[1].id ==~ /^.+_[0-9]+_[0-1]\.[0-9]{2}$/ ?
                        "${it[0].id}_haplotype_${it[1].id.split('_')[1]}" : it[1].id,
                    sequence: it[1].sequence
                ]
            ]
        }
        .map{ ">${it[1].id}\n${it[1].sequence}" }
        .set{ ch_renamed_haplotype }

    consensus_fasta
        .splitFasta(record: [sequence: true], elem: 1)
        .map{ ">${it[0].id}_consensus\n${it[1].sequence}" }
        .set{ ch_renamed_consensus }

    genome_strain           // [accession, strain]
        .join(genome_fasta) // [accession, strain, fasta]
        .map{ it.drop(1) }  // [strain, fasta]
        .map{ [it[0].contains('ROOT') ? 'ROOT' : it[0], it[1]] }
        .splitFasta(record: [sequence: true], elem: 1) // [strain, [sequence]]
        .map{ ">${it[0]}\n${it[1].sequence}" }
        .set{ ch_renamed_genome }

    ch_renamed_haplotype
        .mix(ch_renamed_consensus)
        .mix(ch_renamed_genome)
        .collectFile(name: 'sequences.fasta')
        .map{ [[id: 'collective', single_end: null, strandedness: null], it] }
        .set{ ch_all_sequences }

    MAFFT(ch_all_sequences)
    versions = versions.mix(MAFFT.out.versions)

    RAXMLNG_PARSE(MAFFT.out.fas.map{ it[1] })
    versions = versions.mix(RAXMLNG_PARSE.out.versions)

    RAXMLNG_SEARCH(RAXMLNG_PARSE.out.rba)
    versions = versions.mix(RAXMLNG_SEARCH.out.versions)

    RAXMLNG_BOOTSTRAP(RAXMLNG_PARSE.out.rba)
    versions = versions.mix(RAXMLNG_BOOTSTRAP.out.versions)

    RAXMLNG_SUPPORT(RAXMLNG_SEARCH.out.best_tree, RAXMLNG_BOOTSTRAP.out.bootstraps)
    versions = versions.mix(RAXMLNG_SUPPORT.out.versions)

    tree = RAXMLNG_SUPPORT.out.support

    emit:
    tree
    versions
}
