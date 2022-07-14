include { BLAST_BLASTN } from '../../modules/nf-core/modules/blast/blastn/main.nf'
include { BLAST_MAKEBLASTDB } from '../../modules/nf-core/modules/blast/makeblastdb/main.nf'

workflow CLOSEST_REFERENCE {
    take:
    consensus_fasta
    genome_strain
    genome_fasta

    main:
    versions = Channel.empty()

    // Make a BLAST database out of the strain reference genomes
    BLAST_MAKEBLASTDB(
        genome_fasta
            .map{ it[1] }
            .collectFile(name: 'genomes.fasta', newLine: true)
    )
    versions = versions.mix(BLAST_MAKEBLASTDB.out.versions)

    // BLAST the consensus sequence against all of the reference genomes
    BLAST_BLASTN(consensus_fasta, BLAST_MAKEBLASTDB.out.db.first())
    BLAST_BLASTN.out.txt.map{ [it[0], it[1].readLines()[0]] }.set{ accession }
    versions = versions.mix(BLAST_BLASTN.out.versions)

    // Get the strain name of each sample's closest BLAST hit
    // [meta, strain name]
    accession
        .map{ [it[1], it[0]] }
        .combine(genome_strain, by: 0)
        .map{ [it[1], it[2]] }
        .set{ strain }

    // Get the genome of each sample's closest BLAST hit in fasta format
    accession
        .combine(genome_fasta.map{ [it[1], it[0]] }, by: 1)
        .map{ [it[1], it[2]] }
        .set{ fasta }

    emit:
    accession
    strain
    fasta
    versions
}
