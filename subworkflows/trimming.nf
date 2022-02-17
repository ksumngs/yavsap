#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NANOFILT } from '../modules/local/modules/nanofilt/main.nf'
include { TRIMMOMATIC } from '../modules/local/modules/trimmomatic/main.nf'

workflow TRIMMING {
    take:
    reads

    main:
    if (params.platform == 'illumina') {
        TRIMMOMATIC(reads)
        TRIMMOMATIC.out.fastq.set{ fastq }
        TRIMMOMATIC.out.log.set{ log_out }
    }
    else if (params.platform == 'nanopore') {
        NANOFILT(reads)
        NANOFILT.out.fastq.set{ fastq }
        NANOFILT.out.log.set{ log_out }
    }

    emit:
    fastq
    log_out
}


process read_trimming_ont {
    label 'nanofilt'
    label 'process_low'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), path("${sampleName}_trimmed.fastq.gz"), emit: trimmedreads
    path("${sampleName}.nanofilt.log"), optional: true, emit: report

    script:
    minlenflag = ( params.trim_minlen > 0 )   ? "--length ${params.trim_minlen}"     : ''
    maxlenflag = ( params.trim_maxlen > 0 )   ? "--maxlength ${params.trim_maxlen}"  : ''
    qualflag   = ( params.trim_meanqual > 0 ) ? "--quality ${params.trim_meanqual}"  : ''
    mingcflag  = ( params.trim_mingc > 0 )    ? "--minGC ${params.trim_mingc}"       : ''
    maxgcflag  = ( params.trim_maxgc > 0 )    ? "--maxGC ${params.trim_maxgc}"       : ''
    headflag   = ( params.trim_headcrop > 0 ) ? "--headcrop ${params.trim_headcrop}" : ''
    tailflag   = ( params.trim_tailcrop > 0 ) ? "--tailcrop ${params.trim_tailcrop}" : ''
    optionflags = [
        minlenflag,
        maxlenflag,
        qualflag,
        mingcflag,
        maxgcflag,
        headflag,
        tailflag
    ].join(' ')
    """
    gunzip < ${readsFiles} | \
        NanoFilt --logfile ${sampleName}.nanofilt.log ${optionflags} | \
        gzip > ${sampleName}_trimmed.fastq.gz
    """
}

// Trim Illumina reads
process read_trimming_pe {
    label 'trimmomatic'
    label 'process_medium'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), path("*.fastq.gz"), emit: trimmedreads
    path("${sampleName}.trimmomatic.log"), emit: report

    script:
    // Put together the trimmomatic parameters
    clipflag =
        ( !(params.trim_adapters.getClass() == Boolean || params.trim_adapters.allWhitespace) &&
        params.trim_mismatches > 0 && params.trim_pclip > 0 && params.trim_clip) ?
        "ILLUMINACLIP:/Trimmomatic-0.39/adapters/${params.trim_adapters}:${params.trim_mismatches}:${params.trim_pclip}:${params.trim_clip}" : ''
    winflag =
        ( params.trim_winsize > 0 && params.trim_winqual > 0 ) ? "SLIDINGWINDOW:${params.trim_winsize}:${params.trim_winqual}" : ""
    leadflag =
        ( params.trim_leading > 0 ) ? "LEADING:${params.trim_leading}" : ""
    trailflag =
        ( params.trim_trailing > 0 ) ? "TRAILING:${params.trim_trailing}" : ""
    cropflag =
        ( params.trim_tailcrop > 0 ) ? "CROP:${params.trim_crop}" : ""
    headflag =
        ( params.trim_headcrop > 0 ) ? "HEADCROP:${params.trim_headcrop}" : ""
    minlenflag =
        ( params.trim_minlen > 0 ) ? "MINLEN:${params.trim_minlen}" : ""
    trimsteps = [
        clipflag,
        winflag,
        leadflag,
        trailflag,
        cropflag,
        headflag,
        minlenflag
    ].join(' ')
    """
    trimmomatic PE -threads ${task.cpus} \
        ${readsFiles} \
        ${sampleName}_trimmed_R1.fastq.gz \
        /dev/null \
        ${sampleName}_trimmed_R2.fastq.gz \
        /dev/null \
        ${trimsteps} 2> ${sampleName}.trimmomatic.log
    """
}
