#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENOME_DOWNLOAD } from '../subworkflows/reference.nf'

workflow simulated_reads {
    main:
    GENOME_DOWNLOAD()
    ReferenceGenome = GENOME_DOWNLOAD.out.indexedFasta

    HaplotypeYaml = file("${workflow.projectDir}/test/test.haplotypes.yaml")

    VARIANT_SIMULATE(ReferenceGenome, HaplotypeYaml)
    HaplotypeGenomes = VARIANT_SIMULATE.out

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

process VARIANT_SIMULATE {
    label 'haplink'

    input:
    file(reference)
    file(haplotypeYaml)

    output:
    path "haplotypes.fasta"

    script:
    """
    make-haplotype-fastas ${haplotypeYaml} ${reference[0]} haplotypes.fasta
    """
}

process KRAKEN_READ_SIMULATE {
    label 'kraken'

    input:
    file(genomes)

    output:
    tuple val('simulatedsample'), file("simulatedreads*.fastq.gz")

    script:
    """
    /kraken2-2.1.2/data/simulator.pl --num-frags 4000 \
        --output-format 'simulatedreads_R#.fastq' --read-length 150 \
        --frag-dist-params '500,50' --error_rate 0.01 \
        *.fasta
    gzip *.fastq
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
