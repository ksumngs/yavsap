include { BLAST_BLASTN } from '../modules/nf-core/modules/blast/blastn/main.nf'
include { BLAST_MAKEBLASTDB } from '../modules/local/modules/blast/makeblastdb/main.nf'
include { CAT_CAT } from '../modules/nf-core/modules/cat/cat/main.nf'
include { CUSTOM_ALIGNMENT } from './custom-alignment.nf'
include { EFETCH } from '../modules/local/modules/efetch/main.nf'
include { ESEARCH } from '../modules/local/modules/esearch/main.nf'
include { IVAR_CONSENSUS } from '../modules/nf-core/modules/ivar/consensus/main.nf'
include { SAMTOOLS_BAM2FQ } from '../modules/nf-core/modules/samtools/bam2fq/main.nf'

workflow CLOSEST_REFERENCE {
    take:
    reads
    reference
    genome_list

    main:
    // Transform the TSV genome list into a nf-core sample channel
    Channel
        .fromPath(genome_list)
        .splitCsv(sep: '\t')
        .map{ [
            ['id': it[0], 'single_end': null, 'strandedness': null],
            it[1]
        ] }
        .set{ GenomeTable }

    // Transform the accession number into the search query
    GenomeTable
        .map{ [
            it[0],
            'nucleotide',
            it[1]
        ] }
        .set{ GenomeSearchParameters }

    // Search NCBI for the accession numbers
    ESEARCH(GenomeSearchParameters)

    // Download the matching genomes in fasta format
    ESEARCH.out.xml
        .combine(
            Channel.of(
                ['fasta', null]
            )
        )
        .set{ GenomeFetchParameters }
    EFETCH(GenomeFetchParameters)

    // Create a channel containing every strain's genome in fasta format
    GenomeTable
        .join(EFETCH.out.txt)
        .map{ [it[1], it[2]] }
        .set{ GenomeFastas }

    // Combine all of the genomes together into a single fasta file
    CAT_CAT(
        EFETCH.out.txt
            .map{ it[1] }
            .collect(),
        "genomes.fasta"
    )
    CAT_CAT.out.file_out.set{ genome_fasta }

    // Make a BLAST database out of the strain reference genomes
    BLAST_MAKEBLASTDB(
        genome_fasta
            .map {[
                ['id': 'yavsap-genomes', 'single_end': null, 'strandedness': null],
                it
            ]},
        'nucl'
    )

    // Get the consensus sequence of each sample
    IVAR_CONSENSUS(reads, reference, false)

    // BLAST the consensus sequence against all of the reference genomes
    BLAST_BLASTN(
        IVAR_CONSENSUS.out.fasta,
        BLAST_MAKEBLASTDB.out.db.map{ it[1] }
    )

    BLAST_BLASTN.out.txt
        .map{[
            it[0],
            it[1].readLines()[0]
        ]}
        .set{ accession }

    // Get the strain name of each sample's closest BLAST hit
    accession
        .map{ [it[1], it[0]] }
        .combine(
            GenomeTable
                .map{ [it[1], it[0].id] },
                by: 0
        )
        .map{ [it[1], it[2]] }
        .set{ strain }

    // Get the genome of each sample's closest BLAST hit in fasta format
    accession
        .combine(GenomeFastas.map{ [ it[1], it[0] ] }, by: 1)
        .map{ [ it[1], it[2] ] }
        .set{ fasta }

    // Convert the aligned reads back into fastq format (unalign them?)
    SAMTOOLS_BAM2FQ(reads, true)

    // Align the reads to their new reference genome
    CUSTOM_ALIGNMENT(SAMTOOLS_BAM2FQ.out.reads.join(fasta))
    CUSTOM_ALIGNMENT.out.bam.set{ bam }
    CUSTOM_ALIGNMENT.out.bai.set{ bai }

    emit:
    accession
    strain
    fasta
    bam
    bai
    genome_fasta
}
