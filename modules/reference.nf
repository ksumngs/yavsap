#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Declare what we're going to call our reference genome
ReferenceName = 'JEV'

workflow reference_genome_pull {
    main:
    // Pull and index the reference genome of choice
    reference_genome_pull_fasta | reference_genome_index_samtools
    indexedreference = reference_genome_index_samtools.out

    // Pull and annotate the reference genome of choice
    reference_genome_pull_genbank | reference_genome_annotate
    annotatedreference = reference_genome_annotate.out

    emit:
    indexedreference = indexedreference
    annotatedreference = annotatedreference
}

// Get the reference genome in FASTA format
process reference_genome_pull_fasta {
    cpus 1

    output:
    file '*'

    script:
    """
    efetch -db nucleotide -id ${params.genomeId} -format fasta > reference.fasta
    """
}

// Get the reference genome in GenBank format
process reference_genome_pull_genbank {
    cpus 1

    output:
    file '*'

    script:
    """
    efetch -db nucleotide -id ${params.genomeId} -format gb > reference.gb
    """
}

// Index the reference genome for use with Samtools
process reference_genome_index_samtools {
    cpus 1

    input:
    file(reference)

    output:
    tuple file("${ReferenceName}.fasta"), file("*.fai")

    script:
    """
    # Create a reference genome index
    cp ${reference} ${ReferenceName}.fasta
    samtools faidx ${ReferenceName}.fasta
    """
}

// Process the reference genome's feature table into GFF format
process reference_genome_annotate {
    cpus 1

    input:
    file(reference)

    output:
    file "${ReferenceName}.gff"

    shell:
    '''
    seqret -sequence !{reference} -sformat1 genbank -feature -outseq ref.gff -osformat gff -auto
    head -n $(($(grep -n '##FASTA' ref.gff | cut -d : -f 1) - 1)) ref.gff > !{ReferenceName}.gff
    '''
}
