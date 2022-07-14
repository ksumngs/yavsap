include { EDIRECT_EFETCH } from '../../modules/ksumngs/nf-modules/edirect/efetch/main.nf'
include { EDIRECT_ESEARCH } from '../../modules/ksumngs/nf-modules/edirect/esearch/main.nf'

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

    // Create a list of [strain, [accessions...]]. This has to be a list, and
    // not a channel so that it can be iterated
    def genome_list = []
    genomeFile
        .readLines()
        .each{
            ln ->
                def cells = ln.split('\t')
                def strain = cells[0]
                def accessions = cells.drop(1)
                genome_list << [strain, accessions]
        }

    // Create a list of [accessions...]. This also has to be a list, and is used
    // to sort the sequences within fastas.
    def all_accessions = []
    genome_list
        .each{
            a ->
                a[1].each{
                    b ->
                        all_accessions << b
                }
        }

    // Transform the genome list into a channel
    Channel
        .fromList(genome_list)
        .set{ ch_genome_list }

    // Set the strains output channel
    ch_genome_list
        .map{ [SubworkflowsDownload.accessionTransform(it[1]), it[0]] }
        .set{ strain }
    strain.dump(tag: 'genome_download_strain')

    // Transform the TSV genome list into an edirect query
    ch_genome_list
        .flatMap{ it[1] }
        .toList()
        .map{ it.join(' OR ') }
        .set{ ch_genome_query }
    ch_genome_query.dump(tag: 'genome_download_genome_query')

    // Search NCBI for the accession numbers
    EDIRECT_ESEARCH(ch_genome_query, 'nucleotide')
    versions = versions.mix(EDIRECT_ESEARCH.out.versions)

    // Download the matching genomes in fasta format
    EDIRECT_EFETCH(EDIRECT_ESEARCH.out.xml, 'fasta', '')
    versions = versions.mix(EDIRECT_EFETCH.out.versions)

    /*
        So...
        Another convoluted way to concatenate genome sequences. This is even
        more complicated b/c we have to separate into multiple records. How it
        works:
            1. Start with the fasta file from efetch
            2. Split it into individual records
            3. Sort the records into the same order they were listed in the
                genomes file. Concatenate them into per-strain multi-sequence
                fasta files based on the combined accession numbers
            4. Parse the accession numbers out of the file names
            5. Reformat sequence into accession identifier fasta while removing
                all identifiers and whitespace
            6. Write the results to new files
            7. Extract the accession number out of the file names again
    */
    EDIRECT_EFETCH.out.txt
        .splitFasta(record: [id: true, text: true])
        .collectFile(
            sort: { SubworkflowsDownload.sortSequence(it, all_accessions) }
        ) {
            [
                "${SubworkflowsDownload.accessionKey(it.id, genome_list).replace('|', '-')}.fasta",
                it.text
            ]
        }
        .map{ [it.getBaseName().replaceAll('-', '|'), it] }
        .map{
            [
                it[0],
                ">${it[0]}\n${it[1].text.replaceAll(/>.*\R/, '').replaceAll(/\R/, '')}"
            ]
        }
        .collectFile(){ ["${it[0].replace('|', '-')}.fasta", it[1] ] }
        .map{ [it.getBaseName().replaceAll('-', '|'), it] }
        .set{ fasta }
    fasta.dump(tag: 'genome_fasta')

    emit:
    strain // [accession, strain]
    fasta // [accession, fasta]
    versions
}
