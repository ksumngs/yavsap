#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow trimming {
    take:
    reads

    main:
    if (params.ont) {
        read_trimming_ont(reads)
        trimmedreads = read_trimming_ont.out
    }
    else {
        read_trimming_pe(reads)
        trimmedreads = read_trimming_pe.out
    }

    emit:
    trimmedreads = trimmedreads
}


process read_trimming_ont {
    cpus params.threads

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("${sampleName}_trimmed.fastq.gz")

    script:
    """
    filtlong --min_length ${params.trimMinlen} \
        --keep_percent ${params.trimKeepPercent} \
        --target_bases ${params.trimTargetBases} \
        ${readsFiles} | gzip > ${sampleName}_trimmed.fastq.gz
    """
}

// Trim Illumina reads
process read_trimming_pe {
    cpus params.threads

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("*.fastq.gz")

    script:
    // Put together the trimmomatic parameters
    ILLUMINACLIP = "ILLUMINACLIP:/Trimmomatic-0.39/adapters/${params.trimAdapters}:${params.trimMismatches}:${params.trimPclip}:${params.trimClip}"
    SLIDINGWINDOW = ( params.trimWinsize > 0 && params.trimWinqual > 0 ) ? "SLIDINGWINDOW:${params.trimWinsize}:${params.trimWinqual}" : ""
    LEADING = ( params.trimLeading > 0 ) ? "LEADING:${params.trimLeading}" : ""
    TRAILING = ( params.trimTrailing > 0 ) ? "TRAILING:${params.trimTrailing}" : ""
    CROP = ( params.trimCrop > 0 ) ? "CROP:${params.trimCrop}" : ""
    HEADCROP = ( params.trimHeadcrop > 0 ) ? "HEADCROP:${params.trimHeadcrop}" : ""
    MINLEN = ( params.trimMinlen > 0 ) ? "MINLEN:${params.trimMinlen}" : ""
    trimsteps = ILLUMINACLIP + ' ' + SLIDINGWINDOW + ' ' + LEADING + ' ' + TRAILING + ' ' + CROP + ' ' + HEADCROP + ' ' + MINLEN
    """
    trimmomatic PE -threads ${params.threads} \
        ${readsFiles} \
        ${sampleName}_trimmed_R1.fastq.gz \
        /dev/null \
        ${sampleName}_trimmed_R2.fastq.gz \
        /dev/null \
        ${trimsteps}
    """
}
