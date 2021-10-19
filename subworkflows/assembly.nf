#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow assembly {
    take:
    InputReads
    ReferenceGenome
    GenomeSize

    main:
    if (params.ont) {
        denovo_canu(InputReads, GenomeSize)
        Contigs = denovo_canu.out
    }
    else {
        denovo_spades(InputReads)
        Contigs = denovo_spades.out
    }
    align_to_reference(Contigs, ReferenceGenome)
    AlignedContigs = align_to_reference.out

    emit:
    Contigs
    AlignedContigs
}


// Assemble using Canu
process denovo_canu {
    label 'canu'
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/sequence", mode: "${params.publish_dir_mode}"

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

process denovo_spades {
    label 'spades'
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}/assembly/sequence", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("${sampleName}.contigs.fasta")

    script:
    if (params.spades_mode.getClass() == Boolean || params.spades_mode.allWhitespace) {
        modeflag = ''
    }
    else {
        modeflag = "--${params.spades_mode}"
    }
    """
    spades.py ${modeflag} -o ${sampleName} -1 ${readsFiles[0]} -2 ${readsFiles[1]} -t ${task.cpus}
    cp ${sampleName}/contigs.fasta ./${sampleName}.contigs.fasta
    """

}

// Remap contigs
process align_to_reference {
    label 'minimap'
    publishDir "${params.outdir}/assembly/alignment", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(contigs)
    file(reference)

    output:
    tuple val(sampleName), file("*.{bam,bai}")

    script:
    """
    minimap2 -at ${task.cpus} --MD ${reference[0]} ${contigs} | \
        samtools sort > ${sampleName}.contigs.bam
    samtools index ${sampleName}.contigs.bam
    """
}
