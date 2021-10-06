#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow simulated_reads {
    main:
    JEVGenomes = Channel.from([
        'MT134112.1',
        'M55506.1',
        'L48961.1',
        'MN544779.1',
        'MH753133.1',
        'KR908703.1',
        'NC_014574.1', // Mosquito contamination
        'NC_000845.1'  // Swine contamination
    ])

    download_fasta(JEVGenomes)
    GenomeFastas = download_fasta.out.collect()
    simulate_reads(GenomeFastas)
    OutputReads = simulate_reads.out

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

process simulate_reads {
    label 'kraken'
    label 'process_low'

    input:
    file(genomes)

    output:
    tuple val('simulatedsample'), file("simulatedreads*.fastq.gz")

    script:
    readsoptions = (params.pe) ? \
        "--num-frags 100 --output-format 'simulatedreads_R#.fastq'  --read-length 150  --frag-dist-params '500,50'     --error_rate 0.01" : \
        "--num-frags 347 --output-format 'simulatedreads.fastq'     --read-length 3500 --frag-dist-params '11000,1100' --error_rate 0.0631"
    """
    /kraken2-2.1.2/data/simulator.pl ${readsoptions} *.fasta
    gzip *.fastq
    """
}
