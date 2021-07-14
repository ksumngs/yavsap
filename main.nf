#!/usr/bin/env nextflow

if (params.help) {
    log.info \
    """
NAME
    jev-analysis-pipeline - Automated analysis of Japanese Encephalitis Virus next-generation sequencing data

SYNOPSIS
    nextflow run millironx/jev-analysis-pipeline
        --kraken-db <kraken2 database location>

OPTIONS
    --readsfolder
        The folder containing parired-end Illumina reads in gzipped fastq format. Defaults
        to the current directory

    --threads
        Number of threads to process each sample with. Can't be adjusted on a per-process
        basis. Defaults to 4

    --runname
        A friendly identifier to describe the samples being analyzed. Defaults to
        'jev-analysis'

    --outfolder
        The place where the final anlysis products will be stored. Defaults to runname_out

    --dev
        Run using fewer inputs and faster process options

    --devinputs
        The number of inputs to take in when using --dev

PROCESS-SPECIFIC OPTIONS
Kraken:
    See https://github.com/DerrickWood/kraken2/wiki/Manual for full documentation of
    Kraken 2's available options
    --kraken-db
        Path to Kraken 2 database. REQUIRED
"""
exit 0
}

// Make params persist that need to
RunName = params.runname

// Create an outfolder name if one wasn't provided
if(params.outfolder == "") {
    OutFolder = RunName + "_out"
}
else {
    OutFolder = params.outfolder
}

// Bring in the reads files
Channel
    .fromPath("${params.readsfolder}/*.{fastq,fq}.gz")
    .take( params.dev ? params.devinputs : -1 )
    .map{ file -> tuple(file.simpleName, file) }
    .set{ RawReads }

// Classify reads using Kraken
process kraken {
    cpus params.threads

    input:
    set val(sampleName), file(readsFile) from RawReads

    output:
    tuple sampleName, file("${sampleName}.kraken"), file("${sampleName}.krpt"), file(readsFile) into KrakenFile

    script:
    quickflag = params.dev ? '--quick' : ''
    """
    kraken2 --db ${params.krakenDb} --threads ${params.threads} ${quickflag} \
        --report "${sampleName}.krpt" \
        --output "${sampleName}.kraken" \
        ${readsFile}
    """
}

// Pull the viral reads and any unclassified reads from the original reads
// files for futher downstream processing using KrakenTools
process filterreads {
    cpus 1

    input:
    set val(sampleName), file(krakenFile), file(krakenReport), file(readsFile) from KrakenFile

    output:
    tuple sampleName, file("${sampleName}_filtered.fastq.gz") into FilteredReads

    // Although I haven't seen it documented anywhere, 0 is unclassified reads
    // and 10239 is viral reads
    script:
    """
    extract_kraken_reads.py -k ${krakenFile} \
        -s ${readsFile} \
        -r ${krakenReport} \
        -t 0 10239 --include-children \
        --fastq-output \
        -o ${sampleName}_filtered.fastq
    gzip ${sampleName}_filtered.fastq
    """

}

// Assemble using Canu
process assembly {
    cpus params.threads

    input:
    set val(sampleName), file(readsFile) from FilteredReads

    output:
    tuple val(sampleName), val(assembler), file("${sampleName}.contigs.fasta"), file(readsFile) into ContigsForRemapping

    script:
    assembler = 'metavelvet'
    """
    canu -p ${sampleName} -d out \
        genomeSize=10976 -nanopore ${readsFile}
    cp out/${sampleName}.contigs.fasta .
    """
}

// Get the reference genome
process reference {
    cpus 1
    publishDir OutFolder, mode: 'symlink'

    output:
    file '*' into ReferenceGenome

    script:
    """
    efetch -db nucleotide -id NC_001437 -format fasta > jev.fasta
    """
}

// Remap contigs using BWA
process bwa {
    cpus params.threads

    input:
    set val(sampleName), val(assembler), file(contigs), file(readsFile) from ContigsForRemapping

    output:
    tuple val(sampleName), val(assembler), file(contigs), file("${sampleName}_${assembler}.sam") into RemappedReads

    script:
    """
    cp ${readsFile} read.fastq.gz
    gunzip read.fastq.gz
    bwa index ${contigs}
    bwa aln -t ${params.threads} ${contigs} read.fastq > ${sampleName}.sai
    bwa samse ${contigs} \
        ${sampleName}.sai ${readsFile} > ${sampleName}_${assembler}.sam
    """
}

// Sort and compress the sam files for visualization
process sortsam {
    cpus 1

    input:
    set val(sampleName), val(assembler), file(contigs), file(samfile) from RemappedReads

    output:
    tuple file("${sampleName}_${assembler}.contigs.fasta"), file("${sampleName}_${assembler}.contigs.fasta.fai"), file("${sampleName}_${assembler}.bam"), file("${sampleName}_${assembler}.bam.bai") into Assemblies

    script:
    """
    # Rename the contigs file to a consistent format
    mv ${contigs} ${sampleName}_${assembler}.contigs.fasta

    # Create a contigs indes
    samtools faidx ${sampleName}_${assembler}.contigs.fasta

    # Convert and sort the sam file
    samtools view -S -b ${samfile} > sample.bam
    samtools sort sample.bam -o ${sampleName}_${assembler}.bam

    # Index the sorted bam file
    samtools index ${sampleName}_${assembler}.bam
    """
}

// Create a viewer of all the assembly files
process assemblyview {
    cpus 1

    publishDir OutFolder, mode: 'copy'

    input:
    file '*' from Assemblies.collect()

    output:
    file 'index.html'
    file 'index.js'
    file 'package.json'
    file 'data/*'

    script:
    """
    mkdir data
    mv *.contigs.fasta *.contigs.fasta.fai *.bam *.bam.bai data
    git clone https://github.com/MillironX/igv-bundler.git igv-bundler
    mv igv-bundler/{index.html,index.js,package.json} .
    """
}
