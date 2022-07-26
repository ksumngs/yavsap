include { BLAST_BLASTN } from '../../modules/nf-core/modules/blast/blastn/main.nf'
include { BLAST_MAKEBLASTDB } from '../../modules/nf-core/modules/blast/makeblastdb/main.nf'

workflow CLOSEST_REFERENCE {
    take:
    consensus_fasta
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
    BLAST_BLASTN.out.txt
        .map{
            def accession_num = it[1].readLines()[0]
            def meta = it[0]
            return [meta, accession_num]
        }
        .combine(
            genome_fasta.map{
                def accession_num = it[0]['accession_num']
                def meta = it[0]
                def fasta = it[1]
                return [meta, accession_num, fasta]
            },
            by: 1
        )
        .map{
            def new_meta = it[1].clone()
            new_meta['strain'] = it[2]['id']
            new_meta['blast_accession'] = it[2]['accession_num']
            return [new_meta, it[3]]
        }
        .set{ fasta }
    fasta.dump(tag: 'blast_accession')
    versions = versions.mix(BLAST_BLASTN.out.versions)

    emit:
    fasta
    versions
}
