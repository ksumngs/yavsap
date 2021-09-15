#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow assembly {
    take:
    reads

    main:
    if (params.ont) {
        assembly_ont(reads)
        results = assembly_ont.out
    }
    else {
        assembly_pe(reads)
        results = assembly_pe.out
    }

    emit:
    results
}


// Assemble using Canu
process assembly_ont {
    label 'canu'
    label 'process_medium'

    input:
    tuple val(sampleName), file(readsFile)

    output:
    tuple val(sampleName), file("${sampleName}.contigs.fasta")

    script:
    """
    canu -p ${sampleName} \
        genomeSize=10.5k \
        maxThreads=${task.cpus} \
        maxMemory=${task.memory}
        -nanopore ${readsFile}
    """
}

process assembly_pe {
    label 'spades'
    label 'process_medium'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("contigs.fasta")

    script:
    """
    rnaviralspades.py -o out -1 ${readsFiles[0]} -2 ${readsFiles[1]} -t ${task.cpus}
    cp out/contigs.fasta .
    """

}
