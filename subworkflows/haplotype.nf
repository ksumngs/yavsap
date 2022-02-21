#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PHYLOGENETIC_TREE } from './phylogenetics.nf'
include { CLIQUESNV_ILLUMINA_VC } from '../modules/local/modules/cliquesnv/illumina-vc/main.nf'
include { CLIQUESNV_ILLUMINA } from '../modules/local/modules/cliquesnv/illumina/main.nf'
include { JSON2YAML } from '../modules/local/json2yaml.nf'

workflow HAPLOTYPING {
    take:
    alignments
    references

    main:
    if (params.platform == 'illumina') {
        // Drop the BAM index: CliqueSNV doesn't need it
        alignments
            .map{ it.dropRight(1) }
            .set{ UnindexedAlignments }

        // Do variant calling
        CLIQUESNV_ILLUMINA_VC(UnindexedAlignments)
        CLIQUESNV_ILLUMINA_VC.out.vcf.set{ vcf }

        // Do haplotype calling
        CLIQUESNV_ILLUMINA(UnindexedAlignments)
        CLIQUESNV_ILLUMINA.out.fasta.set{ fasta }

        // Convert haplotyp JSON to YAML
        JSON2YAML(CLIQUESNV_ILLUMINA.out.json)
        JSON2YAML.out.yaml.set{ yaml }
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

    // AllHapSequences = HaplotypeSequences.join(ConsensusSequences)

    // alignment(AllHapSequences, StrainGenomes) | \
    //     PHYLOGENETIC_TREE

    // trees = PHYLOGENETIC_TREE.out

    emit:
    vcf
    yaml
    fasta
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
    mode = (params.platform == 'illumina') ? 'snv-illumina' : 'snv-pacbio'
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
