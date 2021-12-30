#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENOME_DOWNLOAD } from '../subworkflows/reference.nf'

workflow simulated_reads {
    main:
    GENOME_DOWNLOAD()
    ReferenceGenome = GENOME_DOWNLOAD.out.indexedFasta

    HaplotypeYaml = file("${workflow.projectDir}/test/test.haplotypes.yaml")

    VARIANT_SIMULATOR(ReferenceGenome, HaplotypeYaml)
    HaplotypeGenomes = VARIANT_SIMULATOR.out

    if (params.pe) {
        AllGenomes = ReferenceGenome
            .concat(HaplotypeGenomes)
            .collect()
        KRAKEN_READ_SIMULATE(AllGenomes)
        OutputReads = KRAKEN_READ_SIMULATE.out
    }
    else {
        NANOSIM_MODEL_DOWNLOAD()
        NANOSIM_SIMULATE(ReferenceGenome, HaplotypeGenomes, NANOSIM_MODEL_DOWNLOAD.out)
        OutputReads = NANOSIM_SIMULATE.out
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

process VARIANT_SIMULATOR {
    label 'haplink'

    input:
    file(reference)
    file(haplotypeYaml)

    output:
    path "haplotypes.fasta"

    script:
    """
    haplink sequences \
        --haplotypes ${haplotypeYaml} \
        --reference ${reference[0]} \
        --output haplotypes.fasta
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
    /kraken2-2.1.2/data/simulator.pl \
        --num-frags ${depth} \
        --output-format '${prefix}_R#.fastq' \
        --read-length 150 \
        --frag-dist-params '500,50' \
        --error_rate 0.01 \
        ${genomes}
    """
}

process NANOSIM_MODEL_DOWNLOAD {
    label 'run_local'

    output:
    path('metagenome_ERR3152366_Log')

    script:
    """
    wget -qO- \
        https://raw.githubusercontent.com/bcgsc/NanoSim/v3.0.2/pre-trained_models/metagenome_ERR3152366_Log.tar.gz | \
        tar xvz
    """
}

process NANOSIM_SIMULATE {
    label 'nanosim'

    input:
    file(reference)
    file('haplotypes.fasta')
    path(model)

    output:
    tuple val('simulatedsample'), file("simulatedsample.fastq.gz")

    script:
    """
    cp ${reference[0]} REFERENCE.fasta
    cp -v ${workflow.projectDir}/test/{genome_list,abundances,dna_types}.tsv .
    simulator.py metagenome \
        -t ${task.cpus} \
        -c ${model}/training \
        -gl genome_list.tsv \
        -a abundances.tsv \
        -dl dna_types.tsv \
        --seed 42 \
        -b guppy \
        --fastq
    cat *.fastq > simulatedsample.fastq
    gzip simulatedsample.fastq
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
    pbsim --depth ${depth} --hmm_model ${model} ${genome}
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
        cat *_R1.fastq > SIM_R1.fastq
        cat *_R2.fastq > SIM_R2.fastq
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
