#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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

if (!params.ont && !params.pe) {
    log.error "ERROR: Either --ont or --pe must be specified"
    exit 1
}

// Declare what we're going to call our reference genome
ReferenceName = 'JEV'

// Create an output folder name if one wasn't provided
if(params.outfolder == "") {
    OutFolder = params.runname + "_out"
}
else {
    OutFolder = params.outfolder
}

workflow {
    // Pull and index the reference genome of choice
    reference_genome_pull | (reference_genome_index_bowtie & reference_genome_index_samtools)

    // Bring in the reads files
    if (params.ont) {
        raw_reads = Channel
            .fromPath("${params.readsfolder}/*.{fastq,fq}.gz")
            .take( params.dev ? params.devinputs : -1 )
            .map{ file -> tuple(file.simpleName, file) }
    }
    else {
        raw_reads = Channel
            .fromFilePairs("${params.readsfolder}/*{R1,R2,_1,_2}*.{fastq,fq}.gz")
            .take( params.dev ? params.devinputs : -1 )
    }

    // Classify the reads
    raw_reads | read_trimming | read_classification

    // Filter out the non-viral reads
    read_filtering(read_trimming.out, read_classification.out)

    // _de novo_ assemble the viral reads
    assembly(read_filtering.out, reference_genome_pull.out) | \
        contigs_convert_to_fastq

    // Realign contigs to the reference genome
    contigs_realign_to_reference(contigs_convert_to_fastq.out, reference_genome_index_bowtie.out)

    // Realign reads to the reference genome
    reads_realign_to_reference(raw_reads, reference_genome_index_bowtie.out)

    // Make alignments suitable for IGV
    alignment_sort_and_index(contigs_realign_to_reference.out.concat(reads_realign_to_reference.out))

    // Put a pretty bow on everything
    presentation_generator(reference_genome_index_samtools.out, alignment_sort_and_index.out.collect())
}

// Main workflow: will be promoted to ont workflow someday
workflow assembly {
    take:
    reads
    ref

    main:
    if (params.ont) {
        assembly_ont(reads)
        results = assembly_ont.out
    }
    else {
        assembly_pe(reads)
        assembly_pe_improvement(reads, assembly_pe.out, ref)
        results = assembly_pe_improvement.out
    }

    emit:
    results
}

// Get the reference genome
process reference_genome_pull {
    cpus 1

    output:
    file '*'

    script:
    """
    efetch -db nucleotide -id ${params.genomeId} -format fasta > reference.fasta
    """
}

// Index the reference genome
process reference_genome_index_bowtie {
    cpus params.threads

    input:
    file genome

    output:
    file("*.bt2")

    script:
    """
    bowtie2-build --threads ${params.threads} ${genome} ${ReferenceName}
    """
}

process reference_genome_index_samtools {
    cpus 1

    input:
    file(reference)

    output:
    tuple file("${ReferenceName}.fasta"), file("*.fai")

    script:
    """
    # Create a reference genome index
    cp ${reference} ${ReferenceName}.fasta
    samtools faidx ${ReferenceName}.fasta
    """
}

// Trim Illumina reads
process read_trimming {
    cpus params.threads
    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("${sampleName}_trimmed*.fastq.gz")

    script:
    if (params.ont) {
        // Bypass this step for nanopore reads
        """
        cp ${readsFiles} ${sampleName}_trimmed.fastq.gz
        """
    }
    else {
        // Put together the trimmomatic parameters
        ILLUMINACLIP = "ILLUMINACLIP:${params.trimAdapters}:${params.trimMismatches}:${params.trimPclip}:${params.trimClip}"
        SLIDINGWINDOW = ( params.trimWinsize > 0 && params.trimWinqual > 0 ) ? "SLIDINGWINDOW:${params.trimWinsize}:${params.trimWinqual}" : ""
        LEADING = ( params.trimLeading > 0 ) ? "LEADING:${params.trimLeading}" : ""
        TRAILING = ( params.trimTrailing > 0 ) ? "TRAILING:${params.trimTrailing}" : ""
        CROP = ( params.trimCrop > 0 ) ? "CROP:${params.trimCrop}" : ""
        HEADCROP = ( params.trimHeadcrop > 0 ) ? "HEADCROP:${params.trimHeadcrop}" : ""
        MINLEN = ( params.trimMinlen > 0 ) ? "MINLEN:${params.trimMinlen}" : ""
        trimsteps = ILLUMINACLIP + ' ' + SLIDINGWINDOW + ' ' + LEADING + ' ' + TRAILING + ' ' + CROP + ' ' + HEADCROP + ' ' + MINLEN
        """
        trimmomatic PE -threads ${params.threads} \
            ${readsFiles} \
            ${sampleName}_trimmed_R1.fastq.gz \
            /dev/null \
            ${sampleName}_trimmed_R2.fastq.gz \
            /dev/null \
            ${trimsteps}
        """
    }
}

// Classify reads using Kraken
process read_classification {
    cpus params.threads

    input:
    tuple val(sampleName), file(readsFile)

    output:
    tuple file("${sampleName}.kraken"), file("${sampleName}.kreport")

    script:
    quickflag = params.dev ? '--quick' : ''
    pairedflag = params.pe ? '--paired' : ''
    """
    kraken2 --db ${params.krakenDb} --threads ${params.threads} ${quickflag} \
        --report "${sampleName}.kreport" \
        --output "${sampleName}.kraken" \
        ${pairedflag} \
        ${readsFile}
    """
}

// Pull the viral reads and any unclassified reads from the original reads
// files for futher downstream processing using KrakenTools
process read_filtering {
    cpus 1

    input:
    tuple val(sampleName),  file(readsFile)
    tuple file(krakenFile), file(krakenReport)

    output:
    tuple val(sampleName), file("${sampleName}_filtered*.fastq.gz")

    script:
    read2flagin  = (params.pe) ? "-s2 ${readsFile[1]}" : ''
    read1flagout = (params.pe) ? "-o ${sampleName}_filtered_R1.fastq" : " -o ${sampleName}_filtered.fastq"
    read2flagout = (params.pe) ? "-o2 ${sampleName}_filtered_R2.fastq" : ''
    """
    extract_kraken_reads.py -k ${krakenFile} \
        -s ${readsFile[0]} ${read2flagin} \
        -r ${krakenReport} \
        -t ${params.taxIdsToKeep} --include-children \
        --fastq-output \
        ${read1flagout} ${read2flagout}
    gzip -k ${sampleName}_filtered*.fastq
    """
}

// Assemble using Canu
process assembly_ont {
    cpus params.threads

    input:
    tuple val(sampleName), file(readsFile)

    output:
    tuple val(sampleName), file("${sampleName}.contigs.fasta")

    script:
    """
    canu -p ${sampleName} -d out \
        genomeSize=10976\
        maxThreads=${params.threads} \
        stopOnLowCoverage=3 \
        -nanopore ${readsFile}
    cp out/${sampleName}.contigs.fasta .
    """
}

process assembly_pe {
    cpus 1

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    file "contigs.fa"

    script:
    """
    velveth out ${params.kmerLength} -fastq.gz -shortPaired -separate ${readsFiles}
    velvetg out -exp_cov auto -ins_length ${params.velvetInsLength} -min_contig_lgth ${params.velvetMinContigLen}
    cp out/contigs.fa .
    """

}

process assembly_pe_improvement {
    cpus 1

    input:
    tuple val(sampleName), file(readsFiles)
    file(contigs)
    file(reference)

    output:
    tuple val(sampleName), file('scaffolds.scaffolded.gapfilled.length_filtered.sorted.fa')

    script:
    """
    /opt/Bio_AssemblyImprovement-2021.01.04.08.24.00.726/bin/improve_assembly \
        -d -a ${contigs} -f ${readsFiles[0]} -r ${readsFiles[1]} -c ${reference} \
        -s /opt/sspace_basic/SSPACE_Basic.pl \
        -g /usr/local/bin/GapFiller \
        -b /usr/local/bin/abacas.pl
    """
}

// Convert the contigs to fastq with dummy read scores for realignment
process contigs_convert_to_fastq {
    cpus 1

    input:
    tuple val(sampleName), file(contigs)

    output:
    tuple val(sampleName), file("*.fastq.gz")

    script:
    """
    fastx-converter -i ${contigs} -o ${sampleName}.contigs.fastq.gz
    """
}

// Remap contigs using bowtie2
process contigs_realign_to_reference {
    cpus params.threads

    input:
    tuple val(sampleName), file(contigs)
    file(reference)

    output:
    file("${sampleName}.contigs.sam")

    script:
    """
    bowtie2 --threads ${params.threads} -x ${ReferenceName} -U ${contigs} > ${sampleName}.contigs.sam
    """
}

process reads_realign_to_reference {
    cpus params.threads

    input:
    tuple val(sampleName), file(readsFile)
    file(reference)

    output:
    file("${sampleName}.sam")

    script:
    if (params.pe) {
        """
        bowtie2 --threads ${params.threads} -x ${ReferenceName} -1 ${readsFile[0]} -2 ${readsFile[1]} > ${sampleName}.sam
        """
    }
    else {
        """
        bowtie2 --threads ${params.threads} -x ${ReferenceName} -U ${readsFile} > ${sampleName}.sam
        """
    }

}

// Sort and compress the sam files for visualization
process alignment_sort_and_index {
    cpus 1

    input:
    file(samfile)

    output:
    file("*.{bam,bai}")

    script:
    bamfile = samfile.getName().replace('.sam', '.bam')
    """
    # Convert, sort and index the reads file
    samtools view -S -b ${samfile} > sample.bam
    samtools sort sample.bam -o ${bamfile}
    samtools index ${bamfile}

    # Remove intermediate files
    rm sample.bam
    """
}

// Create a viewer of all the assembly files
process presentation_generator {
    cpus 1

    publishDir OutFolder, mode: 'copy'

    input:
    file '*'
    file '*'

    output:
    file 'index.html'
    file 'index.js'
    file 'package.json'
    file 'data/*'

    script:
    """
    mkdir data
    mv *.fasta *.fasta.fai *.bam *.bam.bai data
    git clone https://github.com/MillironX/igv-bundler.git igv-bundler
    mv igv-bundler/{index.html,index.js,package.json} .
    """
}
