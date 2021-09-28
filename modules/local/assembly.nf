#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow assembly {
    take:
    reads
    size

    main:
    if (params.ont) {
        assembly_ont(reads, size)
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
    label 'error_ignore'

    input:
    tuple val(sampleName), file(readsFile)
    val(genomeSize)

    output:
    tuple val(sampleName), file("${sampleName}.contigs.fasta")

    script:
    """
    canu -p ${sampleName} \
        genomeSize=${genomeSize} \
        maxThreads=${task.cpus} \
        useGrid=false \
        correctedErrorRate=${params.canu_corrected_error_rate} \
        minReadLength=${params.canu_min_read_length} \
        minOverlapLength=${params.canu_min_overlap_length} \
        stopOnLowCoverage=${params.canu_stop_on_low_coverage} \
        -nanopore ${readsFile}
    """
}

process assembly_pe {
    label 'spades'
    label 'process_medium'
    label 'error_ignore'

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
