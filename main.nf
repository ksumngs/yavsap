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

include { reference_genome_pull } from './modules/reference.nf'
include { trimming } from './modules/trimming.nf'
include { assembly } from './modules/assembly.nf'

workflow {
    reference_genome_pull()
    IndexedReference = reference_genome_pull.out.indexedreference
    AnnotatedReference = reference_genome_pull.out.annotatedreference

    // Bring in the reads files
    if (params.ont) {
        RawReads = Channel
            .fromPath("${params.readsfolder}/*.{fastq,fq}.gz")
            .take( params.dev ? params.devinputs : -1 )
            .map{ file -> tuple(file.simpleName, file) }
    }
    else {
        RawReads = Channel
            .fromFilePairs("${params.readsfolder}/*{R1,R2,_1,_2}*.{fastq,fq}.gz")
            .take( params.dev ? params.devinputs : -1 )
    }

    RawReads | sample_rename | trimming | read_classification

    KrakenReads = trimming.out.join(read_classification.out)

    // Filter out the non-viral reads
    read_filtering(KrakenReads)
    FilteredReads = read_filtering.out

    // _de novo_ assemble the viral reads
    assembly(FilteredReads)
    Assemblies = assembly.out

    // Realign contigs to the reference genome
    contigs_realign_to_reference(Assemblies, IndexedReference)
    AlignedContigs = contigs_realign_to_reference.out

    // Realign reads to the reference genome
    reads_realign_to_reference(FilteredReads, IndexedReference) | \
        alignments_remove_invalids | \
        alignments_index_samtools
    Alignments = alignments_remove_invalids.out.join(alignments_index_samtools.out)
        .map{n -> [n[0], [n[1], n[2]]] }

    // Call variants
    variants_calling_ivar(Alignments, IndexedReference, AnnotatedReference)
    VariantCalls = variants_calling_ivar.out

    AlignmentStats = Alignments.join(VariantCalls)

    // Get variant stats
    variants_analysis(AlignmentStats, IndexedReference)
    VariantStats = VariantCalls.join(variants_analysis.out)

    // Filter variants
    variants_filter(VariantStats)
    FilteredVariantCalls = Alignments.join(variants_filter.out)

    // Sanity-check those variants
    haplotype_calling_julia(FilteredVariantCalls)

    // Call haplotypes
    haplotype_calling_cliquesnv(Alignments)

    AllAlignments = Alignments.join(AlignedContigs, remainder: true).flatMap{ n -> [n[1], n[2]] }.collect()

    // Put a pretty bow on everything
    presentation_generator(IndexedReference, AllAlignments)
}

process sample_rename {
    cpus 1

    input:
    tuple val(givenName), file(readsFiles)

    output:
    tuple val(sampleName), file("out/*.fastq.gz")

    script:
    sampleName = givenName.split('_')[0]
    if (params.ont) {
        """
        mkdir out
        mv ${readsFiles[0]} out/${sampleName}.fastq.gz
        """
    }
    else {
        """
        mkdir out
        mv ${readsFiles[0]} out/${sampleName}_R1.fastq.gz
        mv ${readsFiles[1]} out/${sampleName}_R2.fastq.gz
        """
    }
}

// Classify reads using Kraken
process read_classification {
    label 'kraken'

    cpus params.threads

    input:
    tuple val(sampleName), file(readsFile)

    output:
    tuple val(sampleName), file("${sampleName}.kraken"), file("${sampleName}.kreport")

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
    label 'krakentools'

    cpus 1

    input:
    tuple val(sampleName), file(readsFile), file(krakenFile), file(krakenReport)

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

// Remap contigs
process contigs_realign_to_reference {
    label 'minimap'

    cpus params.threads

    input:
    tuple val(sampleName), file(contigs)
    file(reference)

    output:
    tuple val(sampleName), file("*.{bam,bai}")

    script:
    """
    minimap2 -at ${params.threads} --MD ${reference[0]} ${contigs} | \
        samtools sort > ${sampleName}.contigs.bam
    samtools index ${sampleName}.contigs.bam
    """
}

process reads_realign_to_reference {
    label 'minimap'

    cpus params.threads

    input:
    tuple val(sampleName), file(readsFile)
    file(reference)

    output:
    tuple val(sampleName), file("${sampleName}_raw.bam")

    script:
    minimapMethod = (params.pe) ? 'sr' : 'map-ont'
    """
    minimap2 -ax ${minimapMethod} -t ${params.threads} --MD ${reference[0]} ${readsFile} | \
        samtools sort > ${sampleName}_raw.bam
    """
}

process alignments_remove_invalids {
    label 'julia'

    cpus 1

    input:
    tuple val(sampleName), file(bamFile)

    output:
    tuple val(sampleName), file("${sampleName}.bam")

    script:
    """
    remove-invalid-reads ${bamFile} ${sampleName}.bam
    """
}

process alignments_index_samtools {
    label 'samtools'

    cpus 1

    input:
    tuple val(sampleName), file(bamFile)

    output:
    tuple val(sampleName), file("*.bai")

    script:
    """
    samtools index ${bamFile}
    """
}

process variants_calling_ivar {
    label 'ivar'

    cpus 1

    input:
    tuple val(sampleName), file(bamfile)
    file(reference)
    file(annotations)

    output:
    tuple val(sampleName), file("*.ivar.tsv")

    script:
    // Crank up the quality metrics (just do it less if we're working with nanopore reads)
    qualFlags = (params.pe) ? '-q30 -t 0.05 -m 1000' : '-q 21 -t 0.05 -m 1500'
    """
    samtools mpileup -aa -A -B -Q 0 --reference ${reference[0]} ${bamfile[0]} | \
        ivar variants -p ${sampleName}.ivar -r ${reference[0]} -g ${annotations} ${qualFlags}
    """
}

// Get stats on the called variants
process variants_analysis {
    label 'bam_readcount'

    cpus 1

    input:
    tuple val(prefix), file(bamfile), file(variantCalls)
    file(reference)

    output:
    tuple val(prefix), file("*.counts.tsv")

    shell:
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
    label 'julia'

    cpus params.threads

    publishDir OutFolder, mode: 'symlink'

    input:
    tuple val(prefix), file(variantCalls), file(variantStats)

    output:
    tuple val(prefix), file("*.filtered.tsv")

    script:
    """
    export JULIA_NUM_THREADS=${params.threads}
    variantfilter ${variantCalls[0]} ${variantStats[0]} ${prefix}.filtered.tsv
    """
}

// At some point, we will need to use long reads to find if mutations are linked within
// a single viral genome. To start, we will look to see if there are reads that contain more
// than one mutation in them as called by ivar
process haplotype_calling_julia {
    label 'julia'

    cpus params.threads

    publishDir OutFolder, mode: 'symlink'

    input:
    tuple val(prefix), file(bamfile), file(variants)

    output:
    file("*.csv")
    file("*.yaml")

    script:
    """
    export JULIA_NUM_THREADS=${params.threads}
    haplotype-finder ${bamfile[0]} ${variants[0]} ${prefix}.haplotypes.csv ${prefix}.haplotypes.yaml
    """
}

process haplotype_calling_cliquesnv {
    label 'cliquesnv'

    cpus params.threads
    memory params.cliquemem

    publishDir OutFolder, mode: 'symlink'

    input:
    tuple val(sampleName), file(bamfile)

    output:
    file "*.{json,fasta}"

    script:
    mode = (params.ont) ? 'snv-pacbio' : 'snv-illumina'
    """
    java -Xmx${params.cliquemem} -jar /opt/CliqueSNV-2.0.2/clique-snv.jar \
        -m ${mode} -threads ${params.threads} -in ${bamfile[0]}
    mv snv_output/* .
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
    cd igv-bundler
    git checkout jev
    cd ..
    mv igv-bundler/{index.html,index.js,package.json} .
    """
}
