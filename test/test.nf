#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENOME_DOWNLOAD } from '../subworkflows/reference.nf'

workflow SIMULATED_READS {
    main:
    // Get the reference genome
    GENOME_DOWNLOAD()
    ReferenceGenome = GENOME_DOWNLOAD.out.indexedFasta

    // Get the haplotype descriptions as a name-file pair channel
    HaplotypeYamls = Channel
        .fromPath("${workflow.projectDir}/test/haplotypes/*.yaml")
        .map{ file -> tuple(file.simpleName, file) }

    // Add muations to the reference genome
    VARIANT_SIMULATOR(HaplotypeYamls, ReferenceGenome)

    // Calculate the requested depth of each haplotype
    HAPLOTYPE_DEPTH(HaplotypeYamls)

    // Create a channel with the mutated genome and the depth
    HaplotypeGenomes = VARIANT_SIMULATOR.out.join(HAPLOTYPE_DEPTH.out)

    // Simulate the reads and return a samplename/fastq tuple
    if (params.platform == 'illumina') {
        HaplotypeGenomes | KRAKEN_READ_SIMULATE
        READ_CONCAT(KRAKEN_READ_SIMULATE.out.collect())
        OutputReads = READ_CONCAT.out
    }
    else {
        PBSIM_MODEL_DOWNLOAD()
        PBSIM_SIMULATE(HaplotypeGenomes, PBSIM_MODEL_DOWNLOAD.out)
        READ_CONCAT(PBSIM_SIMULATE.out.collect())
        OutputReads = READ_CONCAT.out
    }

    emit:
    OutputReads
}

/// summary: |
///   Convert the `frequency` key in a HapLink.jl output file to a valid (integer)
///   read depth
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: The identifier for this sample
///       - name: haplotypes
///         type: file
///         description: |
///           A HapLink.jl file containing a single haplotype definition. If the
///           file contains more than one haplotype, then the process will only
///           return results for the first haplotype.
/// output:
///   - tuple:
///      - type: val(String)
///        description: The identifier as passed through `prefix`
///      - type: env
///        description: The read depth of the haplotype as an integer
process HAPLOTYPE_DEPTH {
    label 'process_low'

    input:
    tuple val(prefix), file(haplotypes)

    output:
    tuple val(prefix), env(DEPTH)

    script:
    """
    depthof ${haplotypes} > DEPTH
    mapfile DEPTH < DEPTH
    export DEPTH
    """
}

/// summary: Simulate the genome of a haplotype using HapLink.jl
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: The unique id for this mutation pattern
///       - name: haplotypeYaml
///         type: file
///         description: YAML file describing the mutations that make up the haplotype
///   - name: reference
///     type: file
///     description: The reference genome being mutated in fasta format
/// output:
///   - tuple:
///       - type: val(String)
///         description: The identifier as passed through `prefix`
///       - type: path
///         description: The mutated sequence in fasta format
process VARIANT_SIMULATOR {
    label 'haplink'

    input:
    tuple val(prefix), file(haplotypeYaml)
    file(reference)

    output:
    tuple val(prefix), path("${prefix}.fasta")

    script:
    """
    haplink sequences \
        --haplotypes ${haplotypeYaml} \
        --reference ${reference[0]} \
        --output ${prefix}.fasta
    """
}

/// summary: Simulate Illumina paired-end reads using Kraken2's perl script
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: The unique id for this mutation pattern
///       - name: genome
///         type: file
///         description: Fasta file containing the sequence to simulate reads of
///       - name: depth
///         type: val(Int)
///         description: The number of reads to simulate
/// output:
///   - type: path
///     description: The fastq files containing the simulated reads
process KRAKEN_READ_SIMULATE {
    label 'kraken'

    input:
    tuple val(prefix), file(genomes), val(depth)

    output:
    path("${prefix}*.fastq")

    script:
    """
    /kraken2-2.1.2/data/simulator.pl \\
        --num-frags ${depth} \\
        --output-format '${prefix}#.fastq' \\
        --read-length 150 \\
        --frag-dist-params '500,50' \\
        --error_rate 0.01 \\
        --random-seed ${params.seed} \\
        ${genomes}
    """
}

/// summary: Download the PBSim model for ONT 9.5 chemistry
/// output:
///   - type: path
///     description: The model file
process PBSIM_MODEL_DOWNLOAD {
    label 'run_local'
    label 'process_low'

    output:
    path('R95.model')

    script:
    """
    wget https://raw.githubusercontent.com/yukiteruono/pbsim2/eaae5b1313e453e5738c591772070ed529b0fad3/data/R95.model
    """
}

/// summary: Simulate ONT long reads using pbsim
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: The unique id for this mutation pattern
///       - name: genome
///         type: file
///         description: Fasta file containing the sequence to simulate reads of
///       - name: depth
///         type: val(Int)
///         description: The number of reads to simulate
///   - name: model
///     type: file
///     description: The HMM model file to simulate with
/// output:
///   - type: path
///     description: The fastq files containing the simulated reads
process PBSIM_SIMULATE {
    label 'pbsim'

    input:
    tuple val(prefix), file(genome), val(depth)
    file(model)

    output:
    path("${prefix}.fastq")

    script:
    """
    pbsim \\
        --depth ${depth} \\
        --hmm_model ${model} \\
        --seed ${params.seed} \\
        ${genome}
    mv sd_0001.fastq ${prefix}.fastq
    """
}

/// summary: |
///   Concatenate and reformat simulated reads into the format generated by a file
///   input read channel
/// input:
///   - name: NA
///     type: file
///     description: All of the reads files
/// output:
///   - tuple:
///       - type: val(String)
///         description: A unique id for this simulated sample
///       - type: path
///         description: The simulated reads files in gzipped fastq format
process READ_CONCAT {
    label 'process_low'

    input:
    file("*")

    output:
    tuple val('SIM'), path("SIM*.fastq.gz")

    script:
    if (params.pe) {
        """
        cat *_1.fastq > SIM_R1.fastq
        cat *_2.fastq > SIM_R2.fastq
        gzip SIM_R?.fastq
        """
    }
    else {
        """
        cat *.fastq > SIM.fastq
        gzip SIM.fastq
        """
    }
}
