#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow reference_genome_pull {
    main:
    // Download the files
    download_fasta()
    download_genbank()
    fastaref   = download_fasta.out
    genbankref = download_genbank.out

    // Measure the size of the genome based on the genbank record
    gensize = genbankref
        .splitText()
        .first()
        .map{ l -> l.replaceAll(/\s+/, ',') }
        .splitCsv()
        .flatten()
        .take(3)
        .last()

    // Get the name of the reference genome
    refname = genbankref
        .splitText()
        .take(2)
        .last()
        .map{ s -> s.replace('DEFINITION  ', '') }
        .map{ s -> s.replace(',', '') }
        .map{ s -> s.replace('.', '') }
        .map{ s -> s.replace(' ', '_') }
        .map{ s -> s.replace('/', '_') }
        .map{ s -> s.trim() }

    // Process the files
    indexing(fastaref, refname)
    indexedreference = indexing.out
    annotation(genbankref, refname)
    annotatedreference = annotation.out

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
    label 'error_backoff'

    output:
    file '*'

    script:
    """
    efetch -db nucleotide -id ${params.genome} -format fasta > reference.fasta
    grep -q '[^[:space:]]' reference.fasta || exit 1
    """
}

// Get the reference genome in GenBank format
process download_genbank {
    label 'edirect'
    label 'run_local'
    label 'process_low'
    label 'error_backoff'

    cpus 1

    output:
    file '*'

    script:
    """
    efetch -db nucleotide -id ${params.genome} -format gb > reference.gb
    grep -q '[^[:space:]]' reference.gb || exit 1
    """
}

// Index the reference genome for use with Samtools
process indexing {
    label 'samtools'
    label 'process_low'
    publishDir "${params.outdir}/reference", mode: "${params.publish_dir_mode}"

    input:
    file(reference)
    val(refname)

    output:
    tuple file("${refname}.fasta"), file("*.fai")

    script:
    """
    # Create a reference genome index
    cp ${reference} ${refname}.fasta
    samtools faidx ${refname}.fasta
    """
}

// Process the reference genome's feature table into GFF format
process annotation {
    label 'seqret'
    label 'process_low'
        publishDir "${params.outdir}/reference", mode: "${params.publish_dir_mode}"

    input:
    file(reference)
    val(refname)

    output:
    file "${refname}.gff"

    shell:
    '''
    seqret -sequence !{reference} -sformat1 genbank -feature -outseq ref.gff -osformat gff -auto
    head -n $(($(grep -n '##FASTA' ref.gff | cut -d : -f 1) - 1)) ref.gff > !{refname}.gff
    '''
}
