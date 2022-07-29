include { EDIRECT_EFETCH } from '../../modules/ksumngs/nf-modules/edirect/efetch/main.nf'
include { EDIRECT_ESEARCH } from '../../modules/ksumngs/nf-modules/edirect/esearch/main.nf'
include { EMBOSS_SEQRET as SEQRET_FASTA } from '../../modules/ksumngs/nf-modules/emboss/seqret/main'

workflow GENOME_DOWNLOAD {
    take:
    genome_tsv

    main:
    versions = Channel.empty()

    // Find the strain genomes list
    genomeFile = file("${genome_tsv}", type: 'file')
    if (!genomeFile.toFile().exists()) {
        genomeFile = file(
            "${workflow.projectDir}/genomes/${genome_tsv}.tsv",
            checkIfExists: true,
            type: 'file'
        )
    }

    def num_strains = 0
    def genome_list = []
    genomeFile
        .readLines()
        .each{
            ln ->
                def cells = ln.split('\t')
                def strain = cells[0]
                def accessions = cells.drop(1)
                num_strains = accessions.length
                accessions.eachWithIndex {
                    accession, index ->
                        genome_list << [
                            [id: strain, strand_num: index, accession_num: accession],
                            accession
                        ]
                }
        }
    Channel.fromList(genome_list).set{ ch_genome_query }
    ch_genome_query.dump(tag: 'genome_query')

    // Search NCBI for the accession numbers
    EDIRECT_ESEARCH(ch_genome_query, 'nucleotide')
    versions = versions.mix(EDIRECT_ESEARCH.out.versions)

    EDIRECT_EFETCH(EDIRECT_ESEARCH.out.xml, 'gb', '')
    versions = versions.mix(EDIRECT_EFETCH.out.versions)

    SEQRET_FASTA(EDIRECT_EFETCH.out.txt, 'fasta')
    versions = versions.mix(SEQRET_FASTA.out.versions)

    /*
        So...
        Another convoluted way to concatenate genome sequences. This is even
        more complicated b/c we have to separate into multiple records. How it
        works:
            1. Start with each fasta file from seqret
            2. Extract their contents
            3. Group and sort sequence contents based on strand order
            4. Recreate identifiers and sequence and create a fasta text string
            5. Put them in files
            6. Re-extract the metadata from the files
    */
    SEQRET_FASTA.out.outseq // [meta, fasta]
        .map{
            def meta_new = it[0].clone()
            meta_new['text'] = it[1].text
            return [
                it[0].id,
                meta_new
            ]
        }
        .groupTuple(by: 0, sort: { it['strand_num'] }, size: num_strains)
        .map{
            def id = it[0]
            def records = it[1]
            def accession_nums = []
            def seq_text = ""

            records.each{
                accession_nums << it['accession_num']
                seq_text += it['text'].replaceAll(/>.*\R/, '').replaceAll(/\R/, '')
            }

            def combo_accession_num = accession_nums.join('|')

            def fasta_text =
                """\
                +>${combo_accession_num} ${id}
                +${seq_text}
                +""".stripMargin('+')

            return [
                [id: id, accession_num: combo_accession_num],
                fasta_text
            ]
        }
        .collectFile(){
            ["${it[0]['accession_num'].replace('|', '-')}.fasta", it[1]]
        }
        .map{
            def accession = (it.text =~ />\S+/).findAll()[0].replace('>', '')
            def id = (it.text =~ / .+/).findAll()[0].trim()
            def file = it

            return [
                [id: id, accession_num: accession],
                file
            ]
        }
        .set{ fasta }
    fasta.dump(tag: 'genome_fasta')

    emit:
    fasta // [[id: name, accession-num: accession], fasta]
    versions
}
