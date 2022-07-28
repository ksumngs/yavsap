#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EDIRECT_EFETCH } from '../../modules/ksumngs/nf-modules/edirect/efetch/main.nf'
include { EDIRECT_ESEARCH } from '../../modules/ksumngs/nf-modules/edirect/esearch/main.nf'
include { EMBOSS_SEQRET as SEQRET_FASTA } from '../../modules/ksumngs/nf-modules/emboss/seqret/main'
include { EMBOSS_SEQRET as SEQRET_GFF } from '../../modules/ksumngs/nf-modules/emboss/seqret/main'
include { GFFCAT } from '../../modules/local/gffcat'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main.nf'

workflow REFERENCE_DOWNLOAD {
    take:
    accession

    main:
    versions = Channel.empty()

    // Convert the accession numbers of the reference genome into metadata
    def accession_list = accession.split(',')
    def accession_plus_meta = []
    accession_list.eachWithIndex {
        accession, index ->
            accession_plus_meta << [
                [id: 'REFERENCE', strand_num: index, accession_num: accession],
                accession
            ]
    }
    Channel.fromList(accession_plus_meta).set{ ch_accession_query }
    def num_strains = accession_list.length
    ch_accession_query.dump(tag: 'accession_query')

    EDIRECT_ESEARCH(ch_accession_query, 'nucleotide')
    EDIRECT_EFETCH(EDIRECT_ESEARCH.out.xml, 'gb', '')

    SEQRET_FASTA(EDIRECT_EFETCH.out.txt, 'fasta')
    SEQRET_GFF(EDIRECT_EFETCH.out.txt, 'gff')

    /*
        So...
        A convoluted way to concatenate genome sequences into a single-sequence
        fasta file. How it works:
            1. Start with the fasta files from efetch
            2. Extract their contents
            3. Combine and sort the sequence contents based on strand order
            4. Remove all identifiers and whitespace from the new fasta record
            5. Write the results into a new fasta format
            6. Remap to an nf-core meta channel
    */
    SEQRET_FASTA.out.outseq // [ meta, fasta ]
        .map{
            def meta_new = it[0].clone()
            meta_new['text'] = it[1].text
            return [
                it[0].id,
                meta_new
            ]
        }
        .groupTuple(by: 0, sort: {it['strand_num']}, size: num_strains)
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

            return [[id: id, accession_num: combo_accession_num], fasta_text]
        }
        .collectFile(name: 'reference.fasta'){ it[1] }
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
    fasta.dump(tag: 'reference')

    fasta
        .map{ it[1] }
        .collectFile(
            name: "reference.fasta",
            storeDir: "${params.outdir}/reference"
        )

    fasta
        .map{ it[1] }
        .collectFile(
            name: "reference.fasta",
            storeDir: "${params.outdir}/report"
        )

    SAMTOOLS_FAIDX(fasta)
    SAMTOOLS_FAIDX.out.fai.set{ fai }

    SEQRET_GFF.out.outseq
        .map{
            def meta_new = it[0].clone()
            return [
                it[0].id,
                [meta_new, it[1]]
            ]
        }
        .groupTuple(
            by: 0,
            sort: { it[0]['strand_num'] },
            size: num_strains
        )
        .map{
            def id = it[0]
            def records = it[1]
            def accession_nums = []
            def gff_files = []

            records.each{
                accession_nums << it[0]['accession_num']
                gff_files << it[1]
            }

            def combo_accession_num = accession_nums.join('|')

            return [[id: id, accession_num: combo_accession_num], gff_files]
        }
        .set{ ch_sorted_gff }
    ch_sorted_gff.dump(tag: 'sorted_gff')

    GFFCAT(ch_sorted_gff)
    GFFCAT.out.gff.set{ gff }

    versions = versions.mix(EDIRECT_ESEARCH.out.versions.first())
    versions = versions.mix(EDIRECT_EFETCH.out.versions.first())
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    versions = versions.mix(SEQRET_FASTA.out.versions.first())
    versions = versions.mix(SEQRET_GFF.out.versions.first())

    emit:
    fasta = fasta.first()
    fai = fai.first()
    gff = gff.first()
    versions
}
