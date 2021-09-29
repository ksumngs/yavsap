#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow haplotyping {
    take:
    Alignments
    Assemblies
    ReferenceGenome

    main:
    calling(Alignments)
    HaplotypeSequences = calling.out.haplotypeSequences.join(Assemblies, remainder: true)

    merge_fastas(HaplotypeSequences, ReferenceGenome) | \
        alignment | \
        phylogenetic_tree
}

process calling {
    label 'cliquesnv'
    label 'process_medium'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(bamfile)

    output:
    tuple val(sampleName), path("${sampleName}.fasta"), emit: haplotypeSequences
    tuple val(sampleName), path("${sampleName}.json"), emit: haplotypeData

    script:
    mode = (params.ont) ? 'snv-pacbio' : 'snv-illumina'
    jmemstring = task.memory.toMega() + 'M'
    """
    java -Xmx${jmemstring} -jar /usr/local/share/cliquesnv/clique-snv.jar \
        -m ${mode} -threads ${task.cpus} -in ${bamfile[0]} -tf 0.01 -fdf extended -rm -log
    mv snv_output/* .
    """
}


process merge_fastas {
    label 'process_low'

    input:
    tuple val(sampleName), file(haplotypes), file(assembly)
    file(reference)

    output:
    tuple val(sampleName), file("${sampleName}.haplotypes.fasta")

    shell:
    '''
    # Keep only the first (most complete) consensus sequence
    if [ "$(grep -c '^>' !{assembly})" -gt 1 ]; then
        head !{assembly} -n $(( $(grep -n '^>' !{assembly} | tail -n +2 | head -n1 | awk '{split($0,a,":"); print a[1]}') - 1 )) > consensus.fasta
    else
        cp !{assembly} consensus.fasta
    fi

    # Label the consensus sequences as such
    sed -i "s/>/>CONSENSUS /g" consensus.fasta

    cat !{reference[0]} consensus.fasta !{haplotypes} > !{sampleName}.haplotypes.fasta
    '''
}

process alignment {
    label 'mafft'
    label 'process_low'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    cpus 1

    input:
    tuple val(sampleName), file(haploReads)

    output:
    tuple val(sampleName), file("${sampleName}.haplotypes.fas")

    shell:
    '''
    mafft --auto !{haploReads} > !{sampleName}.haplotypes.fas
    sed -i "s/ .*$//" !{sampleName}.haplotypes.fas
    '''
}

process phylogenetic_tree {
    label 'raxml'
    label 'error_ignore'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(alignedHaplotypes)

    output:
    tuple val(sampleName), file("${sampleName}.tree")

    script:
    """
    raxml-ng --threads ${task.cpus} --prefix ${sampleName} --all --msa ${alignedHaplotypes} --model --bs-trees 1000 GTR+G
    cp ${sampleName}.raxml.support ${sampleName}.tree
    """
}
