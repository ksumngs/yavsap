#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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
///   -tuple:
///     - type: val(String)
///       description: Sample identifier
///     - type: path
///       description: Multi-alignment file in binary format
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
        --threads ${task.cpus}{auto} \\
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
        --threads ${task.cpus}{auto} \\
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
        --threads ${task.cpus}{auto} \\
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
        --threads ${task.cpus}{auto} \\
        --tree "${tree}" \\
        --bs-trees "${bootstraps}" \\
        --prefix "${prefix}"
    cp "${prefix}.raxml.support" "${prefix}.nwk"
    """
}
