#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { skipping_read } from '../lib/skipping-read.nf'

/// summary: |
///   Take reads from the input folder or a samplesheet and reformat them to be
///   single files with clean names
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: Reads files
workflow READS_INGEST {
    main:
    // Sanity check: interleaved reads cannot be single-end
    if (params.interleaved && !params.paired) {
        log.error "ERROR: --interleaved cannot be specified if --paired is false"
        exit 1
    }

    if (params.samplesheet) {
        SampleList = Channel
            .from(file("${params.samplesheet}"))
            .splitCsv(sep: "\t")
            .filter { !(it[0] ==~ /^#.*/) }

        if (params.paired) {
            RawSamples = SampleList
                .map {
                    [
                        it[0],
                        skipping_read(it.drop(1), 2),
                        skipping_read(it.drop(2), 2)
                    ]
                }
        }
        else {
            RawSamples = SampleList
                .map {
                    [
                        it[0],
                        skipping_read(it.drop(1), 1)
                    ]
                }
        }
    }
    else {
        if (params.paired) {
            RawSamples = Channel
                .fromFilePairs("${params.input}/*{R1,R2,_1,_2}*.{fastq,fq,fastq.gz,fq.gz}")
                .map { [ it[0], it[1][0], it[1][1] ] }
        }
        else {
            RawSamples = Channel
                .fromPath("${params.input}/*.{fastq,fq,fastq.gz,fq.gz}")
                .map{ file -> tuple(file.simpleName, file) }
        }
    }

    if (params.paired) {
        PAIRED_PREPROCESS(RawSamples)
        CleanedSamples = PAIRED_PREPROCESS.out
    }
    else {
        SINGLE_PREPROCESS(RawSamples)
        CleanedSamples = SINGLE_PREPROCESS.out
    }

    emit:
    CleanedSamples
}

/// summary: |
///   Takes an interleaved paired-end read file and splits it into two
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: Sample identifier
///       - name: reads
///         type: file
///         description: The interleaved fastq reads file
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: The forward reads file
///       - type: path
///         description: The reverse reads file
process SEQKIT_SPLIT {
    label 'seqkit'

    input:
    tuple val(prefix), file(reads)

    output:
    tuple val(prefix), path("*.part_001.*"), path("*.part_002.*")

    script:
    """
    seqkit split2 "${reads}" -p2 -O . -f
    """
}

/// summary: |
///   Takes a collection of single-end sequencing reads and converts them into a
///   single file with a safe sample and filename
/// input:
///   - tuple:
///       - name: givenName
///         type: val(String)
///         description: Identifier for this sample
///       - name: reads
///         type: file
///         description: Collection of individual reads files
/// output:
///   - tuple:
///       - type: val(String)
///         description: Identifier for this sample, cleaned of any escape characters
///       - type: path
///         description: The single, complete reads file
process SINGLE_PREPROCESS {
    label 'parallelzip'

    input:
    tuple val(givenName), file(reads)

    output:
    tuple val(sampleName), path("${sampleName}.fastq.gz")

    script:
    sampleName = givenName.split('_')[0]
    """
    parallel -j${task.cpus} gunzip -f ::: ./*.gz
    cat ./*.f*q > "${sampleName}.fastq"
    pigz -p${task.cpus} "${sampleName}.fastq"
    """
}

/// summary: |
///   Takes a collection of paired-end sequencing reads and converts them into a
///   single file with a safe sample and filename
/// input:
///   - tuple:
///       - name: givenName
///         type: val(String)
///         description: Identifier for this sample
///       - name: forwardReads
///         type: file
///         description: Collection of individual forward reads files
///       - name: reverseReads
///         type: file
///         description: Collection of individual reverse reads files
/// output:
///   - tuple:
///       - type: val(String)
///         description: Identifier for this sample, cleaned of any escape characters
///       - type: path
///         description: The combines reads files as paired files
process PAIRED_PREPROCESS {
    label 'parallelzip'

    input:
    tuple val(givenName), file(forwardReads), file(reverseReads)

    output:
    tuple val(sampleName), path("${sampleName}*.fastq.gz")

    script:
    sampleName = givenName.split('_')[0]
    """
    parallel -j${task.cpus} gunzip -f ::: ./*.gz
    cat ./*_*1*.f*q > "${sampleName}_R1.fastq"
    cat ./*_*2*.f*q > "${sampleName}_R2.fastq"
    parallel -j2 pigz -p${task.cpus / 2} ::: "${sampleName}"_R*.fastq
    """
}
