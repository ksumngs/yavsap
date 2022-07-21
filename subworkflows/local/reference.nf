#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EDIRECT_EFETCH } from '../../modules/ksumngs/nf-modules/edirect/efetch/main.nf'
include { EDIRECT_ESEARCH } from '../../modules/ksumngs/nf-modules/edirect/esearch/main.nf'
include { EMBOSS_SEQRET as SEQRET_FASTA } from '../../modules/ksumngs/nf-modules/emboss/seqret/main'
include { EMBOSS_SEQRET as SEQRET_GFF } from '../../modules/ksumngs/nf-modules/emboss/seqret/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main.nf'

workflow REFERENCE_DOWNLOAD {
    take:
    accession

    main:
    versions = Channel.empty()

    def accession_query = accession.replaceAll(',', ' OR ')
    def accession_list = accession.split(',')
    def accession_id = accession.replaceAll(',', '|')

    EDIRECT_ESEARCH(accession_query, 'nucleotide')
    EDIRECT_EFETCH(EDIRECT_ESEARCH.out.xml, 'gb', '')

    EDIRECT_EFETCH.out.txt
        .map{ [ [id: 'genbank'], it ] }
        .set{ ch_efetch_genbank }

    SEQRET_FASTA(ch_efetch_genbank, 'fasta')
    SEQRET_GFF(ch_efetch_genbank, 'gff')

    /*
        So...
        A convoluted way to concatenate genome sequences into a single-sequence
        fasta file. How it works:
            1. Start with the fasta file from efetch
            2. Split it into individual records
            3. Sort the records into the same order they were passed into the
                --genome parameter, and concatenate into a multi-sequence fasta
            4. Remove all identifiers and whitespace from the new fasta file
            5. Add the new identifier to the beginning of the new fasta string
            6. Write the results into a new fasta format
    */
    SEQRET_FASTA.out.outseq // [ meta, fasta ]
        .map{ it[1] }
        .splitFasta(record: [text: true])
        .collectFile(
            name: 'references-sorted.fasta',
            sort: { fa ->
                return SubworkflowsDownload.sortSequence(fa, accession_list.toList())
            }
        ) { it.text }
        .map{ it.text.replaceAll(/>.*\R/, '').replaceAll(/\R/, '') }
        .map{ ">${accession_id}\n${it}" }
        .collectFile(name: 'references-concat.fasta')
        .set{ fasta }
    fasta.dump(tag: 'reference')

    SAMTOOLS_FAIDX(
        fasta.map{ [
            ['id': 'reference'],
            it
        ] }
    )
    SAMTOOLS_FAIDX.out.fai.set{ fai }

    versions = versions.mix(EDIRECT_ESEARCH.out.versions)
    versions = versions.mix(EDIRECT_EFETCH.out.versions)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    fasta = fasta.first()
    fai = fai.first()
    versions
}
