include { BLAST_BLASTN } from '../modules/nf-core/modules/blast/blastn/main.nf'
include { BLAST_MAKEBLASTDB } from '../modules/nf-core/modules/blast/makeblastdb/main.nf'
include { CUSTOM_ALIGNMENT } from './custom-alignment.nf'
include { EDIRECT_EFETCH } from '../modules/ksumngs/nf-modules/edirect/efetch/main.nf'
include { EDIRECT_ESEARCH } from '../modules/ksumngs/nf-modules/edirect/esearch/main.nf'
include { IVAR_CONSENSUS } from '../modules/nf-core/modules/ivar/consensus/main.nf'
include { SAMTOOLS_FASTQ } from '../modules/nf-core/modules/samtools/fastq/main.nf'

workflow CLOSEST_REFERENCE {
    take:
    reads
    reference
    genome_list

    main:
    versions = Channel.empty()

    // Transform the TSV genome list into an edirect query
    genomeQuery = genome_list
        .first()
        .readLines()
        .collect{ it.split('\t')[1] }
        .join(' OR ')

    // Search NCBI for the accession numbers
    EDIRECT_ESEARCH(genomeQuery, 'nucleotide')

    // Download the matching genomes in fasta format
    EDIRECT_EFETCH(EDIRECT_ESEARCH.out.xml, 'fasta', '')
    EDIRECT_EFETCH.out.txt.set{ genome_fasta }

    // Make a BLAST database out of the strain reference genomes
    BLAST_MAKEBLASTDB(genome_fasta)

    // Get the consensus sequence of each sample
    IVAR_CONSENSUS(reads, reference, false)
    IVAR_CONSENSUS.out.fasta.set{ consensus_fasta }

    // BLAST the consensus sequence against all of the reference genomes
    BLAST_BLASTN(
        consensus_fasta,
        BLAST_MAKEBLASTDB.out.db
    )

    BLAST_BLASTN.out.txt
        .map{[
            it[0],
            it[1].readLines()[0]
        ]}
        .set{ accession }

    // Create a channel with strain genome information
    // [accession, strain]
    Channel
        .fromPath(genome_list)
        .splitCsv(sep: '\t')
        .map{ [ it[1], it[0] ] }
        .set{ GenomeTable }

    // Get the strain name of each sample's closest BLAST hit
    // [meta, strain name]
    accession
        .map{ [it[1], it[0]] }
        .combine(GenomeTable, by: 0)
        .map{ [it[1], it[2]] }
        .set{ strain }

    // Create a channel containing every strain's genome in fasta format
    // [accession, fasta]
    genome_fasta
        .splitFasta(file: true)
        .map{ [it.readLines()[0].split(' ')[0].replace('>', ''), it] }
        .set{ GenomeFastas }

    // Get the genome of each sample's closest BLAST hit in fasta format
    accession
        .combine(GenomeFastas.map{ [ it[1], it[0] ] }, by: 1)
        .map{ [ it[1], it[2] ] }
        .set{ fasta }

    // Convert the aligned reads back into fastq format (unalign them?)
    SAMTOOLS_FASTQ(reads)

    // Align the reads to their new reference genome
    CUSTOM_ALIGNMENT(SAMTOOLS_FASTQ.out.fastq.join(fasta))
    CUSTOM_ALIGNMENT.out.bam.set{ bam }
    CUSTOM_ALIGNMENT.out.bai.set{ bai }

    versions = versions.mix(EDIRECT_ESEARCH.out.versions)
    versions = versions.mix(EDIRECT_EFETCH.out.versions)
    versions = versions.mix(BLAST_MAKEBLASTDB.out.versions)
    versions = versions.mix(IVAR_CONSENSUS.out.versions)
    versions = versions.mix(BLAST_BLASTN.out.versions)
    versions = versions.mix(SAMTOOLS_FASTQ.out.versions)
    versions = versions.mix(CUSTOM_ALIGNMENT.out.versions)

    emit:
    accession
    strain
    fasta
    bam
    bai
    genome_fasta
    consensus_fasta
    versions
}
