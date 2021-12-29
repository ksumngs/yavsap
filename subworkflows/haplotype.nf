#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NCBI_DOWNLOAD } from '../modules/ncbi-dl.nf'

workflow haplotyping {
    take:
    Reads
    Alignments
    ReferenceGenome
    GenomeAnnotation

    main:
    GenomePath = "${workflow.projectDir}/genomes/${params.genome_list}.tsv"
    GenomeFile = file(GenomePath)
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

    if (params.pe) {
        calling_pe(RealignedReads)
        HaplotypeSequences = calling_pe.out.haplotypeSequences
    }
    else {
        AlignmentsAndGenomes = RealignedReads.join(BlastGenomes)
        calling_ont(AlignmentsAndGenomes)
        HaplotypeSequences = calling_ont.out.haplotype_fasta
    }

    AllHapSequences = HaplotypeSequences.join(ConsensusSequences)

    alignment(AllHapSequences, StrainGenomes) | \
        phylogenetic_tree

    trees = phylogenetic_tree.out

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

process calling_pe {
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
    java -Xmx${jmemstring} -jar /usr/local/share/cliquesnv/clique-snv.jar \
        -m ${mode} -threads ${task.cpus} -in ${bamfile[0]} -tf 0.01 -fdf extended -rm -log
    mv snv_output/* .
    """
}

process calling_ont {
    label 'haplink'
    label 'process_high'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    input:
    tuple val(prefix), file(bamfile), file(reference)

    output:
    tuple val(prefix), path("variants/${prefix}.vcf"), emit: variants
    tuple val(prefix), path("haplotypes/${prefix}.haplotypes.yaml"), emit: haplotype_yaml
    tuple val(prefix), path("haplotypes/${prefix}.haplotypes.fasta"), emit: haplotype_fasta

    script:
    """
    export JULIA_NUM_THREADS=${task.cpus}
    haplink ${bamfile[0]} \
        --reference ${reference[0]} \
        --variants ${prefix}.vcf \
        --prefix ${prefix} \
        --quality ${params.variant_quality} \
        --frequency ${params.variant_frequency} \
        --position ${params.variant_position} \
        --variant-significance ${params.variant_significance} \
        --variant-depth ${params.variant_depth} \
        --haplotype-significance ${params.haplotype_significance} \
        --haplotype-depth ${params.haplotype_depth}
    rm -rf variants haplotypes
    mkdir variants
    mv ${prefix}.vcf variants
    mkdir haplotypes
    mv ${prefix}.yaml haplotypes/${prefix}.haplotypes.yaml
    mv ${prefix}.fasta haplotypes/${prefix}.haplotypes.fasta
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

process phylogenetic_tree {
    label 'raxml'
    label 'error_ignore'
    publishDir "${params.outdir}/phylogenetics", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(alignedHaplotypes)

    output:
    tuple val(sampleName), file("${sampleName}.nwk")

    script:
    """
    raxml-ng --threads ${task.cpus}{auto} --workers auto \
        --prefix ${sampleName} \
        --all --model GTR+G --bs-trees ${params.phylogenetic_bootstraps} \
        --msa ${alignedHaplotypes}
    cp ${sampleName}.raxml.support ${sampleName}.nwk
    """
}
