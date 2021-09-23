#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

if (params.help) {
    cowsay(
    """\
====================================================================================
                                JEV Analysis Pipeline
====================================================================================

    jev-analysis-pipeline - Intra-sample viral population analysis targeted at
    Japanese Encephalitis

    Usage:

        nextflow run millironx/jev-analysis-pipeline

    Options:

        --input             Relative or absolute path to directory containing
                            gzipped fastq files
                                type: path, default: .

        --platform          Type of reads to process. Options are 'illumina' and
                            'nanopore'
                                type: string, default: none

        --genome            NCBI accession number of the reference genome to align
                            reads against
                                type: string, default: 'NC_001437.1'

        --kraken2_db        Kraken2-compatible database for classifying reads
                                type: path, default: none

        --keep_taxid        Space-separated list of NCBI taxids to keep and process
                            after classifying
                                type: string, default: '0 10239'

        --outdir            Directory in which to place results
                                type: path, default: ./results

        --help              Print this message and exit

    For more information on usage and parameters, visit the website at
        https://github.com/MillironX/jev-analysis-pipeline
"""
)
exit 0
}

if (!params.ont && !params.pe) {
    log.error "ERROR: --platform <illumina,nanopore> must be specified"
    exit 1
}

// Declare what we're going to call our reference genome
ReferenceName = 'JEV'

include { reference_genome_pull } from './modules/local/reference.nf'
include { trimming }              from './modules/local/trimming.nf'
include { assembly }              from './modules/local/assembly.nf'

cowsay(
"""\
====================================================================================
                                JEV Analysis Pipeline
====================================================================================

Input folder:           ${params.input}
Sequencing platform:    ${params.platform}
    Illumina?:          ${params.pe}
    Nanopore?:          ${params.ont}
Reference genome:       ${params.genome}
Kraken2 Database:       ${params.kraken2_db}
Taxonomic Ids:          '${params.keep_taxid}'
Output folder           ${params.outdir}
Diagnostics folder:     ${params.tracedir}
"""
)

workflow {
    reference_genome_pull()
    IndexedReference = reference_genome_pull.out.indexedreference
    AnnotatedReference = reference_genome_pull.out.annotatedreference

    // Bring in the reads files
    if (params.ont) {
        RawReads = Channel
            .fromPath("${params.input}/*.{fastq,fq}.gz")
            .map{ file -> tuple(file.simpleName, file) }
    }
    else {
        RawReads = Channel
            .fromFilePairs("${params.input}/*{R1,R2,_1,_2}*.{fastq,fq}.gz")
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
    reads_realign_to_reference(FilteredReads, IndexedReference)
    Alignments = reads_realign_to_reference.out

    /*
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

    Haplotypes = haplotype_calling_julia.out.haplotype_yaml
    ConsensusHaplotypes = Haplotypes.join(Assemblies)
    haplotype_conversion_fasta(ConsensusHaplotypes, IndexedReference) | \
        haplotype_alignment | \
        haplotype_phylogenetic_tree
    */

    AllAlignments = Alignments.join(AlignedContigs, remainder: true).flatMap{ n -> [n[1], n[2]] }.collect()

    haplotype_calling_cliquesnv(Alignments)
    HaplotypeSequences = haplotype_calling_cliquesnv.out.haplotypeSequences.join(Assemblies, remainder: true)

    haplotype_merge_fastas(HaplotypeSequences, IndexedReference) | \
        haplotype_alignment | \
        haplotype_phylogenetic_tree

    // Put a pretty bow on everything
    presentation_generator(IndexedReference, AllAlignments)
}

process sample_rename {
    label 'process_low'

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
    label 'process_high_memory'

    input:
    tuple val(sampleName), file(readsFile)

    output:
    tuple val(sampleName), file("${sampleName}.kraken"), file("${sampleName}.kreport")

    script:
    pairedflag = params.pe ? '--paired' : ''
    """
    kraken2 --db ${params.kraken2_db} --threads ${task.cpus} \
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
    label 'process_low'

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
        -t ${params.keep_taxid} --include-children \
        --fastq-output \
        ${read1flagout} ${read2flagout}
    gzip -k ${sampleName}_filtered*.fastq
    """
}

// Remap contigs
process contigs_realign_to_reference {
    label 'minimap'

    input:
    tuple val(sampleName), file(contigs)
    file(reference)

    output:
    tuple val(sampleName), file("*.{bam,bai}")

    script:
    """
    minimap2 -at ${task.cpus} --MD ${reference[0]} ${contigs} | \
        samtools sort > ${sampleName}.contigs.bam
    samtools index ${sampleName}.contigs.bam
    """
}

process reads_realign_to_reference {
    label 'minimap'

    input:
    tuple val(sampleName), file(readsFile)
    file(reference)

    output:
    tuple val(sampleName), file("${sampleName}.{bam,bam.bai}")

    script:
    minimapMethod = (params.pe) ? 'sr' : 'map-ont'
    """
    minimap2 -ax ${minimapMethod} -t ${task.cpus} --MD ${reference[0]} ${readsFile} | \
        samtools sort > ${sampleName}.bam
    samtools index ${sampleName}.bam
    """
}

process variants_calling_ivar {
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
process variants_analysis {
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
process variants_filter {
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
process haplotype_calling_julia {
    label 'julia'
    label 'error_retry'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    input:
    tuple val(prefix), file(bamfile), file(variants)

    output:
    tuple val(prefix), path("${prefix}.haplotypes.csv"), emit: linkage_stats
    tuple val(prefix), path("${prefix}.haplotypes.yaml"), emit: haplotype_yaml

    script:
    """
    export JULIA_NUM_THREADS=${task.cpus}
    haplotype-finder ${bamfile[0]} ${variants[0]} ${prefix}.haplotypes.csv ${prefix}.haplotypes.yaml
    """
}

process haplotype_conversion_fasta {
    label 'julia'
    label 'error_retry'
    label 'process_low'

    input:
    tuple val(sampleName), file(haplotypeYaml), file(assembly)
    file(reference)

    output:
    tuple val(sampleName), file("${sampleName}.haplotypes.fasta")

    shell:
    '''
    # Keep only the first (most complete) consensus sequence
    cp !{assembly} consensus.fasta

    # Label the consensus sequences as such
    sed -i "s/>/>CONSENSUS /g" consensus.fasta

    make-haplotype-fastas !{haplotypeYaml} !{reference[0]} !{sampleName}.fasta

    '''
}

process haplotype_alignment {
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

process haplotype_phylogenetic_tree {
    label 'raxml'
    label 'error_ignore'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    input:
    tuple val(sampleName), file(alignedHaplotypes)

    output:
    tuple val(sampleName), file("${sampleName}.raxml.bestTree")

    script:
    """
    raxml-ng --threads ${task.cpus} --prefix ${sampleName} --all --msa ${alignedHaplotypes} --model GTR+G
    """
}

process haplotype_calling_cliquesnv {
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

process haplotype_merge_fastas {
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

// Create a viewer of all the assembly files
process presentation_generator {
    label 'process_low'
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

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

def cowsay(message) {
messagelines = message.split('\n')
nlines = messagelines.length
linelength = 0
messagelines.each{ l -> if ( l.length() > linelength ) { linelength = l.length() } }
paddinglength = linelength + 2

if ( nlines == 1 ) {
    balloon =
""" ${"_"*paddinglength}
< ${message} >
 ${"-"*paddinglength}"""
}
else {
balloon =
""" ${"_"*paddinglength}
/ ${messagelines[0].padRight(linelength)} \\"""
for (int i=1;i<(nlines-1);i++) {
balloon =
"""${balloon}
| ${messagelines[i].padRight(linelength)} |"""
}
balloon =
"""${balloon}
\\ ${messagelines[nlines-1].padRight(linelength)} /
 ${"-"*paddinglength}"""
}

cow =
"""        \\   ^__^
         \\  (oo)\\_______
            (__)\\       )\\/\\
                ||----w |
                ||     ||"""

cowput =
"""${balloon}
${cow}
"""

log.info cowput
}
