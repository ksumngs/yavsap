#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow trimming {
    take:
    reads

    main:
    if (params.ont) {
        read_trimming_ont(reads)
        samplename = read_trimming_ont.out.samplename
        trimmedreads = read_trimming_ont.out.trimmedreads
    }
    else {
        read_trimming_pe(reads)
        samplename = read_trimming_pe.out.samplename
        trimmedreads = read_trimming_pe.out.trimmedreads
    }

    emit:
    samplename = samplename
    trimmedreads = trimmedreads
}


process read_trimming_ont {
    cpus params.threads

    input:
    tuple val(fileName), file(readsFiles)

    output:
    val(sampleName), emit: samplename
    path "*.fastq.gz", emit: trimmedreads

    script:
    sampleName = fileName.split('_')[0]
    """
    gunzip -c ${readsFiles} | \
        NanoFilt -l 100 -q 10 | \
        gzip > ${sampleName}_trimmed.fastq.gz
    """

}

// Trim Illumina reads
process read_trimming_pe {
    cpus params.threads
    input:
    tuple val(fileName), file(readsFiles)

    output:
    val(sampleName), emit: samplename
    path "*.fastq.gz", emit: trimmedreads

    script:
    sampleName = fileName.split('_')[0]
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
