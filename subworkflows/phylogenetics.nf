#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CAT_CAT } from '../modules/nf-core/modules/cat/cat/main.nf'
include { CAT_FASTQ } from '../modules/nf-core/modules/cat/fastq/main.nf'
include { MAFFT } from '../modules/local/modules/mafft/main.nf'
include { RAXMLNG_BOOTSTRAP } from '../modules/local/modules/raxmlng/bootstrap/main.nf'
include { RAXMLNG_PARSE } from '../modules/local/modules/raxmlng/parse/main.nf'
include { RAXMLNG_SEARCH } from '../modules/local/modules/raxmlng/search/main.nf'
include { RAXMLNG_SUPPORT } from '../modules/local/modules/raxmlng/support/main.nf'
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

    RAXMLNG_PARSE(MAFFT.out.fas)

    RAXMLNG_SEARCH(RAXMLNG_PARSE.out.rba)
    RAXMLNG_BOOTSTRAP(RAXMLNG_PARSE.out.rba)

    RAXMLNG_SUPPORT(RAXMLNG_SEARCH.out.best_tree.join(RAXMLNG_BOOTSTRAP.out.bootstraps))

    RAXMLNG_SUPPORT.out.support.set{ tree }

    emit:
    tree
}

/// summary: |
///   Convert a plain-text multi-alignment to RAxML-NG's binary alignment format
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
///         description: Multi-alignment file in binary format
process RAXML_PARSE {
    label 'raxml'
    label 'error_ignore'
    label 'process_low'

    input:
    tuple val(prefix), file(alignment)

    output:
    tuple val(prefix), path("${prefix}.raxml.rba")

    script:
    """
    raxml-ng \\
        --parse \\
        --threads auto{${task.cpus}} \\
        --msa "${alignment}" \\
        --model GTR+G \\
        --prefix "${prefix}"
    """
}

/// summary: Infer the best phylogenetic tree of an alignment
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: Sample identifier
///       - name: alignment
///         type: file
///         description: Multi-alignment file to infer tree from
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: Best inferred tree
process RAXML_SEARCH {
    label 'raxml'
    label 'process_high'

    input:
    tuple val(prefix), file(alignment)

    output:
    tuple val(prefix), path("${prefix}.raxml.bestTree")

    script:
    """
    raxml-ng \\
        --threads auto{${task.cpus}} \\
        --workers auto \\
        --msa "${alignment}" \\
        --model GTR+G \\
        --prefix "${prefix}" \\
        --seed ${params.seed}
    """
}

/// summary: Perform phylogenetic tree bootstrapping on an alignment
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: Sample identifier
///       - name: alignment
///         type: file
///         description: Multi-alignment file to perform bootstrapping on
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: Bootstrap iterations of the trees for this alignment
process RAXML_BOOTSTRAP {
    label 'raxml'
    label 'process_high'

    input:
    tuple val(prefix), file(alignment)

    output:
    tuple val(prefix), path("${prefix}.raxml.bootstraps")

    script:
    """
    raxml-ng \\
        --bootstrap \\
        --threads auto{${task.cpus}} \\
        --workers auto \\
        --msa "${alignment}" \\
        --model GTR+G \\
        --prefix "${prefix}" \\
        --bs-trees ${params.phylogenetic_bootstraps} \\
        --bs-cutoff ${params.phylogenetic_bootstrap_cutoff} \\
        --seed ${params.seed}
    """
}

/// summary: Infer the best phylogenetic tree of an alignment
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: Sample identifier
///       - name: tree
///         type: file
///         description: Best inferred tree
///       - name: bootstraps
///         type: path
///         description: Bootstrap iterations of the trees for this alignment
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: Annotated support tree
process RAXML_SUPPORT {
    label 'raxml'
    label 'process_low'
    publishDir "${params.outdir}/phylogenetics", mode: "${params.publish_dir_mode}"

    input:
    tuple val(prefix), file(tree), file(bootstraps)

    output:
    tuple val(prefix), path("${prefix}.nwk")

    script:
    """
    raxml-ng \\
        --support \\
        --threads auto{${task.cpus}} \\
        --tree "${tree}" \\
        --bs-trees "${bootstraps}" \\
        --prefix "${prefix}"
    cp "${prefix}.raxml.support" "${prefix}.nwk"
    """
}
