#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CLIQUESNV_ILLUMINA } from '../modules/local/modules/cliquesnv/illumina/main.nf'
include { CLIQUESNV_ILLUMINA_VC } from '../modules/local/modules/cliquesnv/illumina-vc/main.nf'
include { HAPLINK_HAPLOTYPES } from '../modules/local/modules/haplink/haplotypes/main.nf'
include { HAPLINK_SEQUENCES } from '../modules/local/modules/haplink/sequences/main.nf'
include { HAPLINK_VARIANTS } from '../modules/local/modules/haplink/variants/main.nf'
include { JSON2YAML } from '../modules/local/json2yaml.nf'
include { PHYLOGENETIC_TREE } from './phylogenetics.nf'

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
        HAPLINK_VARIANTS(alignments.join(references))
        HAPLINK_VARIANTS.out.vcf.set{ vcf }

        HAPLINK_HAPLOTYPES(
            alignments
                .map{ it.dropRight(1) }
                .join(HAPLINK_VARIANTS.out.vcf)
        )
        HAPLINK_HAPLOTYPES.out.yaml.set{ yaml }

        HAPLINK_SEQUENCES(
            HAPLINK_HAPLOTYPES.out.yaml
                .join(references)
        )
        HAPLINK_SEQUENCES.out.fasta.set{ fasta }
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
