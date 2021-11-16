#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { reference_genome_pull } from '../subworkflows/reference.nf'

workflow simulated_reads {
    main:
    reference_genome_pull()
    ReferenceGenome = reference_genome_pull.out.indexedreference

    ContaminantGenomes = Channel.from([
        'NC_014574.1',      // Mosquito contamination
        'NC_000845.1'       // Swine contamination
    ])

    download_fasta(ContaminantGenomes)

    simulate_mutations(ReferenceGenome)

    GenomeFastas = simulate_mutations.out
        .splitFasta( file: true )
        .concat(download_fasta.out)
        .collect()

    if (params.pe) {
        simulate_pe_reads(GenomeFastas)
        OutputReads = simulate_pe_reads.out
    }
    else {
        simulate_ont_reads(GenomeFastas)
        OutputReads = simulate_ont_reads.out
    }

    emit:
    OutputReads
}

process download_fasta {
    label 'edirect'
    label 'process_low'
    label 'error_backoff'
    label 'run_local'

    input:
    val(accessionNumber)

    output:
    file("${accessionNumber}.fasta")

    script:
    """
    efetch -db nucleotide -id ${accessionNumber} -format fasta > \
        ${accessionNumber}.fasta
    grep -q '[^[:space:]]' ${accessionNumber}.fasta || exit 1
    """
}

process simulate_mutations {
    label 'julia'
    label 'error_retry'

    input:
    file(reference)

    output:
    path "haplotypes.fasta"

    script:
    """
    cp ${workflow.projectDir}/test/test.haplotypes.yaml .
    make-haplotype-fastas test.haplotypes.yaml ${reference[0]} haplotypes.fasta
    """
}

process simulate_pe_reads {
    label 'kraken'
    label 'process_low'

    input:
    file(genomes)

    output:
    tuple val('simulatedsample'), file("simulatedreads*.fastq.gz")

    script:
    """
    /kraken2-2.1.2/data/simulator.pl --num-frags 100 \
        --output-format 'simulatedreads_R#.fastq' --read-length 150 \
        --frag-dist-params '500,50' --error_rate 0.01 \
        *.fasta
    gzip *.fastq
    """
}

process simulate_ont_reads {
    label 'nanosim'
    label 'run_local'

    input:
    file(genomes)

    output:
    tuple val('simulatedsample'), file("simulatedsample.fastq.gz")

    script:
    """
    wget -qO- https://raw.githubusercontent.com/bcgsc/NanoSim/v3.0.2/pre-trained_models/metagenome_ERR3152366_Log.tar.gz | tar xvz
    mv -v metagenome_ERR3152366_Log/* .
    cp -v ${workflow.projectDir}/test/{genome_list,abundances,dna_types}.tsv .
    simulator.py metagenome -gl genome_list.tsv -a abundances.tsv -dl dna_types.tsv --seed 42 -b guppy --fastq -t ${task.cpus}
    cat *.fastq > simulatedsample.fastq
    gzip simulatedsample.fastq
    """
}
