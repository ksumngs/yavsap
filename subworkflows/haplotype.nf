#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NCBI_DOWNLOAD } from '../modules/ncbi-dl.nf'
include { PHYLOGENETIC_TREE } from './phylogenetics.nf'

workflow haplotyping {
    take:
    Reads
    Alignments
    ReferenceGenome
    GenomeAnnotation

    main:
    GenomePath = params.genome_list
    GenomeFile = file(GenomePath)
    if (!GenomeFile.toFile().exists()) {
        GenomePath = "${workflow.projectDir}/genomes/${params.genome_list}*"
        GenomeFile = file(GenomePath, checkIfExists: true)
    }
    GenomeList = Channel
        .fromPath(GenomePath)
        .splitCsv(sep: '\t')
        .map{ n -> "${n[1]}"}

    NCBI_DOWNLOAD(GenomeList)
    GenomeFastas = NCBI_DOWNLOAD.out.fasta.collect()

    GENOME_LABELING(GenomeFile, GenomeFastas)
    AccessionGenomes = GENOME_LABELING.out.accessionGenomes
    StrainGenomes = GENOME_LABELING.out.strainGenomes

    blast_db(AccessionGenomes)
    BlastDb = blast_db.out

    consensus(Alignments)
    ConsensusSequences = consensus.out

    blast_consensus(ConsensusSequences, BlastDb)
    BlastHits = blast_consensus.out

    AccessionGenomesSequences = AccessionGenomes
        .splitFasta(record: [id: true, seqString: true])
        .map{ f -> [f.id, ">${f.id}\n${f.seqString}"] }

    BlastGenomes = BlastHits
        .map { h -> [h[1], h[0]] }
        .combine(AccessionGenomesSequences, by: 0)
        .map{ g -> g[1..2] }

    ReadsAndGenomes = Reads.join(BlastGenomes)

    realign_to_new_reference(ReadsAndGenomes)

    RealignedReads = realign_to_new_reference.out.alignment

    if (params.platform == 'illumina') {
        CLIQUESNV_VARIANTS(RealignedReads)
        CLIQUESNV_HAPLOTYPES(RealignedReads)
        HaplotypeSequences = CLIQUESNV_HAPLOTYPES.out.haplotypeSequences
    }
    else {
        AlignmentsAndGenomes = RealignedReads.join(BlastGenomes)
        HAPLINK_VARIANTS(AlignmentsAndGenomes)

        AlignmentsAndVariants = RealignedReads.join(HAPLINK_VARIANTS.out)
        HAPLINK_HAPLOTYPES(AlignmentsAndVariants)

        HaplotypesAndGenomes = HAPLINK_HAPLOTYPES.out.join(BlastGenomes)
        HAPLINK_FASTA(HaplotypesAndGenomes)

        HaplotypeSequences = HAPLINK_FASTA.out
    }

    AllHapSequences = HaplotypeSequences.join(ConsensusSequences)

    alignment(AllHapSequences, StrainGenomes) | \
        PHYLOGENETIC_TREE

    trees = PHYLOGENETIC_TREE.out

    emit:
    trees
}

process GENOME_LABELING {
    label 'process_low'

    input:
    file(genomeList)
    path(genomeFastas)

    output:
    path('accession_genomes.fasta'), emit: accessionGenomes
    path('strain_genomes.fasta'), emit: strainGenomes

    script:
    """
    labelstrains \
        ${genomeList} \
        accession_genomes.fasta \
        strain_genomes.fasta \
        ${genomeFastas}
    """
}

process blast_db {
    label 'blast'

    input:
    path(genomes)

    output:
    tuple val(dbname), path("${dbname}.fasta*")

    script:
    dbname = "YAVSAP_${workflow.sessionId}"
    """
    cp ${genomes} ${dbname}.fasta
    makeblastdb -in ${dbname}.fasta -title ${dbname} -dbtype nucl
    """
}

process consensus {
    label 'ivar'

    input:
    tuple val(sampleName), path(bamfile)

    output:
    tuple val(sampleName), path("${sampleName}.consensus.fasta")

    script:
    """
    samtools mpileup -aa -A -d0 -Q 0 ${bamfile[0]} | \
        ivar consensus \
        -p ${sampleName} \
        -q ${params.variant_quality} \
        -t ${params.variant_frequency} \
        -m ${params.variant_depth}
    mv ${sampleName}.fa ${sampleName}.consensus.fasta
    """
}

process blast_consensus {
    label 'blast'

    input:
    tuple val(sampleName), path(consensusSequence)
    tuple val(blastDbName), path(blastdb)

    output:
    tuple val(sampleName), env(TOPBLASTHIT)

    shell:
    '''
    TOPBLASTHIT=$(blastn -query !{consensusSequence} \
        -db !{blastDbName}.fasta \
        -num_alignments 1 \
        -outfmt "6 saccver" \
        -num_threads !{task.cpus} | head -n1)
    '''
}

process realign_to_new_reference {
    label 'minimap'
    publishDir "${params.outdir}/alignment", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), path(reads), file(referenceGenome)

    output:
    tuple val(sampleName), path("${sampleName}.bam{,.bai}"), emit: alignment
    path("${sampleName}_REFERENCE.fasta{,.fai}"), emit: genome

    script:
    minimapMethod = (params.pe) ? 'sr' : 'map-ont'
    """
    cp ${referenceGenome} ${sampleName}_REFERENCE.fasta
    samtools faidx ${sampleName}_REFERENCE.fasta
    minimap2 -ax ${minimapMethod} \
        -t ${task.cpus} \
        --MD ${sampleName}_REFERENCE.fasta \
        ${reads} | \
        samtools sort > ${sampleName}.bam
    samtools index ${sampleName}.bam
    """
}

/// summary: Call variants on Illumina reads using CliqueSNV
/// input:
///   - tuple:
///       - name: sampleName
///         type: val(String)
///         description: Unique identifier for this sample
///       - name: bamfile
///         type: file
///         description: Alignment file to call variants from
/// output:
///   - tuple:
///       - type: val(String)
///         description: The sample identifier passed through `sampleName`
///       - type: path
///         description: Variant calls in VCF format
process CLIQUESNV_VARIANTS {
    label 'cliquesnv'
    label 'process_high'
    publishDir "${params.outdir}/variants", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(bamfile)

    output:
    tuple val(sampleName), path("${sampleName}.vcf")

    script:
    jmemstring = task.memory.toMega() + 'M'
    """
    java -Xmx${jmemstring} -jar /usr/local/share/cliquesnv/clique-snv.jar \\
        -m snv-illumina-vc \\
        -in ${bamfile[0]} \\
        -t ${params.haplotype_depth} \\
        -tf ${params.haplotype_frequency} \\
        -log \\
        -threads ${task.cpus} \\
        -outDir .
    """
}

/// summary: Call haplotypes for Illumina reads using CliqueSNV
/// input:
///   - tuple:
///       - name: sampleName
///         type: val(String)
///         description: Unique identifier for this sample
///       - name: bamfile
///         type: file
///         description: Alignment file to call variants from
/// output:
///   - name: haplotypeSequences
///     tuple:
///       - type: val(String)
///         description: The sample identifier passed through `sampleName`
///       - type: path
///         description: The sequences of the found haplotypes in FASTA format
///   - name: haplotypeData
///     tuple:
///       - type: val(String)
///         description: The sample identifier passed through `sampleName`
///       - type: path
///         description: Descriptions of the found haplotypes in JSON format
process CLIQUESNV_HAPLOTYPES {
    label 'cliquesnv'
    label 'process_medium'
    publishDir "${params.outdir}/haplotypes", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(bamfile)

    output:
    tuple val(sampleName), path("${sampleName}.fasta"), emit: haplotypeSequences
    tuple val(sampleName), path("${sampleName}.json"), emit: haplotypeData

    script:
    mode = (params.ont) ? 'snv-pacbio' : 'snv-illumina'
    jmemstring = task.memory.toMega() + 'M'
    """
    java -Xmx${jmemstring} -jar /usr/local/share/cliquesnv/clique-snv.jar \\
        -m snv-illumina \\
        -in ${bamfile[0]} \\
        -t ${params.haplotype_depth} \\
        -tf ${params.haplotype_frequency} \\
        -log \\
        -threads ${task.cpus} \\
        -outDir . \\
        -fdf extended
    """
}

/// summary: Call variants for Oxford Nanopore reads using HapLink.jl
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: The identifier for this sample as used in the output filename
///       - name: bamfile
///         type: file
///         description: File for the alignment to call variants from
///       - name: reference
///         type: file
///         description: Reference genome to call variants from in fasta format
/// output:
///   - tuple:
///       - type: val(String)
///         description: The identifier as passed though `prefix`
///       - type: path
///         description: The variants found in the reads in VCF format
process HAPLINK_VARIANTS {
    label 'haplink'
    publishDir "${params.outdir}/variants", mode: "${params.publish_dir_mode}"

    input:
    tuple val(prefix), file(bamfile), file(reference)

    output:
    tuple val(prefix), path("${prefix}.vcf")

    script:
    """
    haplink variants \
        --bam ${bamfile[0]} \
        --reference ${reference} \
        --output ${prefix}.vcf \
        --quality ${params.variant_quality} \
        --frequency ${params.variant_frequency} \
        --position ${params.variant_position} \
        --significance ${params.variant_significance} \
        --depth ${params.variant_depth} \
        --julia-args -t${task.cpus}
    """
}

/// summary: Call haplotypes from long reads using HapLink.jl
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: The identifier for this sample as used in the output filename
///       - name: bamfile
///         type: file
///         description: File for the alignment to call haplotypes from
///       - name: variants
///         type: file
///         description: File containing variants to consider haplotypes make of
/// output:
///   - tuple:
///       - type: val(String)
///         description: The identifier as passed though `prefix`
///       - type: path
///         description: The found haplotypes described in YAML format
process HAPLINK_HAPLOTYPES {
    label 'haplink'
    label 'process_high'
    publishDir "${params.outdir}/haplotypes", mode: "${params.publish_dir_mode}"

    input:
    tuple val(prefix), file(bamfile), file(variants)

    output:
    tuple val(prefix), path("${prefix}.haplotypes.yaml")

    script:
    """
    haplink haplotypes \\
        --bam "${bamfile[0]}" \\
        --variants "${variants}" \\
        --output "${prefix}.haplotypes.yaml" \\
        --significance ${params.haplotype_significance} \\
        --depth ${params.haplotype_depth} \\
        --method ${params.haplotype_method} \\
        --overlap-min ${params.haplotype_overlap_min} \\
        --overlap-max ${params.haplotype_overlap_max} \\
        --iterations ${params.haplotype_iterations} \\
        --seed ${params.seed} \\
        --julia-args -t${task.cpus}
    """
}

/// summary: |
///   Convert a YAML describing haplotypes into a fasta file with the sequences
///   of those haplotypes
/// input:
///   - tuple:
///       - name: prefix
///         type: val(String)
///         description: The identifier for this sample as used in the output filename
///       - name: haplotypes
///         type: file
///         description: The YAML file describing the SNPs in each haplotype
///       - name: reference
///         type: file
///         description: Reference genome to mutate sequences in to create the haplotype sequences
/// output:
///   - tuple:
///       - type: val(String)
///         description: The identifier as passed though `prefix`
///       - type: path
///         description: The mutated sequences in fasta format
process HAPLINK_FASTA {
    label 'haplink'
    label 'process_low'
    publishDir "${params.outdir}/haplotypes", mode: "${params.publish_dir_mode}"

    input:
    tuple val(prefix), file(haplotypes), file(reference)

    output:
    tuple val(prefix), path("${prefix}.haplotypes.fasta")

    script:
    """
    haplink sequences \
        --haplotypes ${haplotypes} \
        --reference ${reference} \
        --output ${prefix}.haplotypes.fasta
    """
}

process alignment {
    label 'mafft'
    label 'process_medium'
    publishDir "${params.outdir}/multi_alignment", mode: "${params.publish_dir_mode}"

    cpus 1

    input:
    tuple val(sampleName), file(haploReads), file(consensus)
    path(referenceStrains)

    output:
    tuple val(sampleName), file("${sampleName}.haplotypes.fas")

    script:
    """
    cat ${consensus} ${haploReads} ${referenceStrains} > ${sampleName}.mafft.fasta
    mafft --thread ${task.cpus} --auto \
        ${sampleName}.mafft.fasta > ${sampleName}.haplotypes.fas
    sed -i "s/ .*\$//" ${sampleName}.haplotypes.fas
    """
}
