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

MORE INFO
    https://github.com/MillironX/jev-analysis-pipeline
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
    reference_genome_pull_fasta | (reference_genome_index_bowtie & reference_genome_index_samtools)

    // Pull and annotate the reference genome of choice
    reference_genome_pull_genbank | reference_genome_annotate

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

    // Realign reads to the reference genome
    reads_realign_to_reference(read_filtering.out, reference_genome_index_samtools.out)

    /*
    // _de novo_ assemble the viral reads
    assembly(read_filtering.out, reference_genome_pull_fasta.out) | \
        contigs_convert_to_fastq

    // Realign contigs to the reference genome
    contigs_realign_to_reference(contigs_convert_to_fastq.out, reference_genome_index_bowtie.out)

    // Make alignments suitable for IGV
    alignment_sort_and_index(reads_realign_to_reference.out.concat(contigs_realign_to_reference.out))
    */

    // Call variants
    variants_calling_ivar(reads_realign_to_reference.out, reference_genome_index_samtools.out, reference_genome_annotate.out)

    // Get variant stats
    variants_analysis(variants_calling_ivar.out, reads_realign_to_reference.out, reference_genome_index_samtools.out)

    // Filter variants
    variants_filter(variants_calling_ivar.out, variants_analysis.out)

    // Sanity-check those variants
    multimutation_search(reads_realign_to_reference.out, variants_filter.out)

    // Put a pretty bow on everything
    presentation_generator(reference_genome_index_samtools.out, reads_realign_to_reference.out.collect())
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
        results = assembly_pe.out
    }

    emit:
    results
}

// Get the reference genome in FASTA format
process reference_genome_pull_fasta {
    cpus 1

    output:
    file '*'

    script:
    """
    efetch -db nucleotide -id ${params.genomeId} -format fasta > reference.fasta
    """
}

// Get the reference genome in GenBank format
process reference_genome_pull_genbank {
    cpus 1

    output:
    file '*'

    script:
    """
    efetch -db nucleotide -id ${params.genomeId} -format gb > reference.gb
    """
}

// Index the reference genome for use with Bowtie2
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

// Index the reference genome for use with Samtools
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

// Process the reference genome's feature table into GFF format
process reference_genome_annotate {
    cpus 1

    input:
    file(reference)

    output:
    file "${ReferenceName}.gff"

    shell:
    '''
    seqret -sequence !{reference} -sformat1 genbank -feature -outseq ref.gff -osformat gff -auto
    head -n $(($(grep -n '##FASTA' ref.gff | cut -d : -f 1) - 1)) ref.gff > !{ReferenceName}.gff
    '''
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
        ILLUMINACLIP = "ILLUMINACLIP:/Trimmomatic-0.39/adapters/${params.trimAdapters}:${params.trimMismatches}:${params.trimPclip}:${params.trimClip}"
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
    cpus params.threads

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("contigs.fasta")

    script:
    """
    rnaviralspades.py -o out -1 ${readsFiles[0]} -2 ${readsFiles[1]} -t ${params.threads}
    cp out/contigs.fasta .
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
    seqtk seq -F '~' ${contigs} | gzip > ${sampleName}.contigs.fastq.gz
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
    file("*.{bam,bai}")

    script:
    minimapMethod = (params.pe) ? 'sr' : 'map-ont'
    """
    minimap2 -ax ${minimapMethod} -t ${params.threads} --MD ${reference[0]} ${readsFile} | \
        samtools sort > ${sampleName}.bam
    samtools index ${sampleName}.bam
    """
}

process variants_calling_ivar {
    cpus 1

    publishDir OutFolder, mode: 'symlink'

    input:
    file(bamfile)
    file(reference)
    file(annotations)

    output:
    file("*.ivar.tsv")

    script:
    // We have to refer to the first file in each of the inputs b/c they are tuples
    // containing the required index files as well
    prefix = bamfile[0].getName().replace('.bam', '')

    // Crank up the quality metrics (just do it less if we're working with nanopore reads)
    qualFlags = (params.pe) ? '-q30 -t 0.05 -m 1000' : '-q 21 -t 0.05 -m 1500'
    """
    samtools mpileup -aa -A -B -Q 0 --reference ${reference[0]} ${bamfile[0]} | \
        ivar variants -p ${prefix}.ivar -r ${reference[0]} -g ${annotations} ${qualFlags}
    """
}

// Get stats on the called variants
process variants_analysis {
    cpus 1

    input:
    file(variantCalls)
    file(bamfile)
    file(reference)

    output:
    file("*.counts.tsv")

    shell:
    prefix = bamfile[0].getName().replace('.bam', '')
    '''
    touch !{prefix}.counts.tsv
    while read -r LINE; do
        echo "$LINE" | while IFS=$'\\t' read -r -a CELLS; do
            REGION="${CELLS[0]}"
            POS="${CELLS[1]}"
            if [[ $POS != "POS" ]]; then
                bam-readcount -f !{reference[0]} !{bamfile[0]} "${REGION}:${POS}-${POS}" >> \
                    !{prefix}.counts.tsv
            fi
        done
    done < !{variantCalls}
    '''
}

// More strictly filter the variants based on strand bias and read position
process variants_filter {
    cpus 1

    input:
    file(variantCalls)
    file(variantStats)

    output:
    file("*.filtered.tsv")

    script:
    prefix = variantCalls[0].getName().replace('.tsv', '')
    """
    variantfilter ${variantCalls[0]} ${variantStats[0]} ${prefix}.filtered.tsv
    """
}

// At some point, we will need to use long reads to find if mutations are linked within
// a single viral genome. To start, we will look to see if there are reads that contain more
// than one mutation in them as called by ivar
process multimutation_search {
    cpus 1

    publishDir OutFolder, mode: 'symlink'

    input:
    file(bamfile)
    file(variants)

    output:
    file("*.varreport")
    file("*.csv")

    script:
    prefix = bamfile[0].getName().replace('.bam', '')
    """
    find-variant-reads ${bamfile[0]} ${variants[0]} ${prefix}.matrix.csv > ${prefix}.varreport
    """
}

// Create a viewer of all the assembly files
process presentation_generator {
    cpus 1

    publishDir OutFolder, mode: 'symlink'

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
