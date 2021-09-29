#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Declare what we're going to call our reference genome
ReferenceName = 'JEV'

workflow reference_genome_pull {
    main:
    // Pull and index the reference genome of choice
    download_fasta | indexing
    indexedreference = indexing.out

    // Pull and annotate the reference genome of choice
    download_genbank()
    genbankref = download_genbank.out
    annotation(genbankref)
    annotatedreference = annotation.out

    // Measure the size of the genome based on the genbank record
    gensize = genbankref
        .splitText()
        .first()
        .map{ l -> l.replaceAll(/\s+/, ',') }
        .splitCsv()
        .flatten()
        .take(3)
        .last()

    emit:
    indexedreference = indexedreference
    annotatedreference = annotatedreference
    genomesize = gensize
}

// Get the reference genome in FASTA format
process download_fasta {
    label 'edirect'
    label 'run_local'
    label 'process_low'

    output:
    file '*'

    script:
    """
    efetch -db nucleotide -id ${params.genome} -format fasta > reference.fasta
    """
}

// Get the reference genome in GenBank format
process download_genbank {
    label 'edirect'
    label 'run_local'
    label 'process_low'

    cpus 1

    output:
    file '*'

    script:
    """
    efetch -db nucleotide -id ${params.genome} -format gb > reference.gb
    """
}

// Index the reference genome for use with Samtools
process indexing {
    label 'samtools'
    label 'process_low'

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
process annotation {
    label 'seqret'
    label 'process_low'

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
