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
    cpus params.threads

    input:
    tuple val(sampleName), file(readsFile)

    output:
    tuple val(sampleName), file("${sampleName}.contigs.fasta")

    script:
    """
    canu -p ${sampleName} -d out \
        genomeSize=10976\
        maxThreads=${params.threads} \
        stopOnLowCoverage=3 \
        -nanopore ${readsFile}
    cp out/${sampleName}.contigs.fasta .
    """
}

process assembly_pe {
    cpus params.threads

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("contigs.fasta")

    script:
    """
    rnaviralspades.py -o out -1 ${readsFiles[0]} -2 ${readsFiles[1]} -t ${params.threads}
    cp out/contigs.fasta .
    """

}
