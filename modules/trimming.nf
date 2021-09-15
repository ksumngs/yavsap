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
    label 'filtlong'
    label 'process_low'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("${sampleName}_trimmed.fastq.gz")

    script:
    """
    filtlong --min_length ${params.trim_minlen} \
        --keep_percent ${params.trim_keep_percent} \
        --target_bases ${params.trim_target_bases} \
        ${readsFiles} | gzip > ${sampleName}_trimmed.fastq.gz
    """
}

// Trim Illumina reads
process read_trimming_pe {
    label 'trimmomatic'
    label 'process_medium'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("*.fastq.gz")

    script:
    // Put together the trimmomatic parameters
    ILLUMINACLIP = "ILLUMINACLIP:/Trimmomatic-0.39/adapters/${params.trim_adapters}:${params.trim_mismatches}:${params.trim_pclip}:${params.trim_clip}"
    SLIDINGWINDOW = ( params.trim_winsize > 0 && params.trim_winqual > 0 ) ? "SLIDINGWINDOW:${params.trim_winsize}:${params.trim_winqual}" : ""
    LEADING = ( params.trim_leading > 0 ) ? "LEADING:${params.trim_leading}" : ""
    TRAILING = ( params.trim_trailing > 0 ) ? "TRAILING:${params.trim_trailing}" : ""
    CROP = ( params.trim_crop > 0 ) ? "CROP:${params.trim_crop}" : ""
    HEADCROP = ( params.trim_headcrop > 0 ) ? "HEADCROP:${params.trim_headcrop}" : ""
    MINLEN = ( params.trim_minlen > 0 ) ? "MINLEN:${params.trim_minlen}" : ""
    trimsteps = ILLUMINACLIP + ' ' + SLIDINGWINDOW + ' ' + LEADING + ' ' + TRAILING + ' ' + CROP + ' ' + HEADCROP + ' ' + MINLEN
    """
    trimmomatic PE -threads ${task.cpus} \
        ${readsFiles} \
        ${sampleName}_trimmed_R1.fastq.gz \
        /dev/null \
        ${sampleName}_trimmed_R2.fastq.gz \
        /dev/null \
        ${trimsteps}
    """
}
