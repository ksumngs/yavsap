#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow haplotyping {
    take:
    Alignments
    Assemblies
    ReferenceGenome
    GenomeAnnotation

    main:
    if (params.pe) {
        calling(Alignments)
        HaplotypeSequences = calling_pe.out.haplotypeSequences.join(Assemblies, remainder: true)
    }
    else {
        variant_calling(Alignments, ReferenceGenome, GenomeAnnotation)
        VariantCalls = variant_calling.out

        AlignmentStats = Alignments.join(VariantCalls)

        // Get variant stats
        alignment_analysis(AlignmentStats, ReferenceGenome)
        VariantStats = VariantCalls.join(alignment_analysis.out)

        // Filter variants
        variant_filter(VariantStats)
        FilteredVariantCalls = Alignments.join(variant_filter.out)

        // Sanity-check those variants
        calling_ont(FilteredVariantCalls, ReferenceGenome)

        HaplotypeSequences = calling_ont.out.haplotype_fasta.join(Assemblies, remainder: true)
    }

    merge_fastas(HaplotypeSequences, ReferenceGenome) | \
        alignment | \
        phylogenetic_tree

    trees = phylogenetic_tree.out

    emit:
    trees
}

process calling_pe {
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

process variant_calling {
    label 'ivar'
    label 'process_low'

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
process alignment_analysis {
    label 'bam_readcount'
    label 'process_low'

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
process variant_filter {
    label 'julia'
    label 'error_retry'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    input:
    tuple val(prefix), file(variantCalls), file(variantStats)

    output:
    tuple val(prefix), file("*.filtered.tsv")

    script:
    """
    export JULIA_NUM_THREADS=${task.cpus}
    variantfilter ${variantCalls[0]} ${variantStats[0]} ${prefix}.filtered.tsv
    """
}

// At some point, we will need to use long reads to find if mutations are linked within
// a single viral genome. To start, we will look to see if there are reads that contain more
// than one mutation in them as called by ivar
process calling_ont {
    label 'julia'
    label 'error_retry'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    input:
    tuple val(prefix), file(bamfile), file(variants)
    file(reference)

    output:
    tuple val(prefix), path("${prefix}.haplotypes.yaml"), emit: haplotype_yaml
    tuple val(prefix), path("${prefix}.haplotypes.fasta"), emit: haplotype_fasta

    script:
    """
    export JULIA_NUM_THREADS=${task.cpus}
    haplotype-finder -p ${params.haplotype_significance} -m ${params.haplotype_minimum} \
        ${bamfile[0]} ${variants[0]} ${prefix}.haplotypes.yaml
    make-haplotype-fastas ${prefix}.haplotypes.yaml ${reference[0]} ${prefix}.haplotypes.fasta
    """
}

process merge_fastas {
    label 'process_low'

    input:
    tuple val(sampleName), file(haplotypes), file(assembly)
    file(reference)

    output:
    tuple val(sampleName), file("${sampleName}.population.fasta")

    script:
    """
    # Create an editable copy of the assembly
    cp ${assembly} assembly.fasta

    # Ensure the records have unique names by prepending the sample name
    # This work for contigs generated by SPAdes and Canu
    sed -i "s/>/>${sampleName}./g" assembly.fasta

    # cat them togther
    cat ${reference[0]} assembly.fasta ${haplotypes} > ${sampleName}.population.fasta
    """
}

process alignment {
    label 'mafft'
    label 'process_medium'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    cpus 1

    input:
    tuple val(sampleName), file(haploReads)

    output:
    tuple val(sampleName), file("${sampleName}.haplotypes.fas")

    shell:
    '''
    mafft --thread !{task.cpus} --auto \
        !{haploReads} > !{sampleName}.haplotypes.fas
    sed -i "s/ .*$//" !{sampleName}.haplotypes.fas
    '''
}

process phylogenetic_tree {
    label 'raxml'
    label 'error_ignore'
    publishDir "${params.outdir}/data", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(alignedHaplotypes)

    output:
    tuple val(sampleName), file("${sampleName}.tree")

    script:
    """
    raxml-ng --threads ${task.cpus} --prefix ${sampleName} --all --msa ${alignedHaplotypes} --model GTR+G --bs-trees 1000
    cp ${sampleName}.raxml.support ${sampleName}.tree
    """
}
