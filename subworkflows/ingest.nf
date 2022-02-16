#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { skipping_read } from '../lib/skipping-read.nf'
include { CAT_FASTQ } from '../modules/nf-core/modules/cat/fastq/main.nf'
include { SEQKIT_SPLIT2 } from '../modules/nf-core/modules/seqkit/split2/main.nf'

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
    // First sanity check: --input must exist
    if (!file(params.input).exists()) {
        log.error "ERROR: file or directory '${params.input}' does not exist!"
        exit 1
    }

    // Sanity check: interleaved reads cannot be single-end
    if (params.interleaved && !params.paired) {
        log.error "ERROR: --interleaved cannot be specified if --paired is false"
        exit 1
    }

    if (file(params.input).isFile()) {
        // --input represents a samplesheet

        // Parse the samplesheet to nf-core sample channel
        Channel
            .of(file("${params.input}"))
            .splitCsv(sep: "\t")
            .filter { !(it[0] ==~ /^#.*/) }
            .map {
                [
                    ['id': it[0], 'single_end': !(params.paired && !params.interleaved), 'strandedness': null],
                    skipping_read(it.drop(1), 1)
                ]
            }
            .set { SampleList }

        // Calculate how many reads can be in a sample channel before it needs to be
        // `cat`ed together
        def maxSamples = (params.paired && !params.interleaved) ? 2 : 1

        // Separate the cat-requiring reads from the regular reads
        SampleList
            .branch {
                meta, fastq ->
                    single: fastq.size() == maxSamples
                        return [ meta, fastq.flatten() ]
                    multiple: fastq.size() > maxSamples
                        return [ meta, fastq.flatten() ]
            }
            .set { CountedSamples }

        // Concatenate the reads together
        CAT_FASTQ(CountedSamples.multiple)
            .reads
            .mix(CountedSamples.single)
            .set { InterleavedSamples }
    }
    else if (file(params.input).isDirectory()) {
        // --input represents a directory of reads
        if (params.paired && !params.interleaved) {
            // Paired reads can be directly sent to interleaved jail
            Channel
                .fromFilePairs("${params.input}/*{R1,R2,_1,_2}*.{fastq,fq,fastq.gz,fq.gz}")
                .map { [
                        [ 'id': it[0].split('_')[0], 'single_end': false, 'strandedness': null ],
                        it[1]
                        ]
                    }
                .set { InterleavedSamples }
        }
        else {
            // Reformat single-end/interleaved reads
            Channel
                .fromPath("${params.input}/*.{fastq,fq,fastq.gz,fq.gz}")
                .map{ file -> tuple(file.simpleName, file) }
                .map{ [
                    [ 'id': it[0].split('_')[0], 'single_end': true, 'strandedness': null ],
                    it[1]
                ] }
                .set { InterleavedSamples }
        }
    }

    if (params.interleaved) {
        // Deinterleave any interleaved reads
        SEQKIT_SPLIT2(InterleavedSamples)
        SEQKIT_SPLIT2.out.reads
            .map { [
                ['id': it[0]['id'], 'single_end': false, 'strandedness': null ],
                it[1]
            ] }
            .set { sample_info }
    }
    else {
        // Transfer reads out of interleaved jail and output them
        InterleavedSamples.set { sample_info }
    }

    emit:
    sample_info
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
