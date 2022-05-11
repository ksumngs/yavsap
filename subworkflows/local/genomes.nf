include { EDIRECT_EFETCH } from '../../modules/ksumngs/nf-modules/edirect/efetch/main.nf'
include { EDIRECT_ESEARCH } from '../../modules/ksumngs/nf-modules/edirect/esearch/main.nf'

workflow GENOME_DOWNLOAD {
    take:
    genome_list

    main:
    versions = Channel.empty()

    // Find the strain genomes list
    genomeFile = file("${genome_list}", type: 'file')
    if (!genomeFile.toFile().exists()) {
        genomeFile = file(
            "${workflow.projectDir}/genomes/${genome_list}.tsv",
            checkIfExists: true,
            type: 'file'
        )
    }

    // Transform the genome list into a channel
    Channel
        .fromPath(genomeFile)
        .splitCsv(sep: '\t')
        .map{ [it[1], it[0]] }
        .set{ strain }

    // Transform the TSV genome list into an edirect query
    genomeQuery = genomeFile
        .readLines()
        .collect{ it.split('\t')[1] }
        .join(' OR ')

    // Search NCBI for the accession numbers
    EDIRECT_ESEARCH(genomeQuery, 'nucleotide')
    versions = versions.mix(EDIRECT_ESEARCH.out.versions)

    // Download the matching genomes in fasta format
    EDIRECT_EFETCH(EDIRECT_ESEARCH.out.xml, 'fasta', '')
    EDIRECT_EFETCH
        .out
        .txt
        .splitFasta(file: true)
        .map{ [it.readLines()[0].split(' ')[0].replace('>', ''), it] }
        .set{ fasta }
    versions = versions.mix(EDIRECT_EFETCH.out.versions)

    emit:
    strain // [accession, strain]
    fasta // [accession, fasta]
    versions
}
