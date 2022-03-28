
include { HAPLOTYPE_YAML2TSV } from '../modules/local/haplotype-yaml2tsv.nf'

workflow PRESENTATION {
    take:
    bam
    bai
    reference_fasta
    reference_fai
    consensus_fasta
    accession
    strain
    haplotype_yaml
    haplotype_fasta

    main:
    versions = Channel.empty()

    HAPLOTYPE_YAML2TSV(haplotype_yaml)
    HAPLOTYPE_YAML2TSV.out.tsv
        .map{ it[1] }
        .collectFile(name: 'haplotypes.tsv')
        .set{ haplotype_tsv }

    versions = versions.mix(HAPLOTYPE_YAML2TSV.out.versions)

    consensus_fasta
        .splitFasta(record: [seqString: true])
        .map{ [ it[0], "${it[0].id}_consensus", null, it[1].seqString ] }
        .set{ consensus_table }

    bam
        .join(bai)
        .join(accession)
        .join(strain)
        .map{ [
            ['id': it[0].id, 'single_end': null, 'strandedness': null],
            it[1],
            it[2],
            it[3],
            it[4]
        ] }
        .set{ strain_table }

    haplotype_tsv
        .splitCsv(sep: '\t')
        .map{ [
            ['id': it[0], 'single_end': null, 'strandedness': null],
            it[1],
            it[2]
        ] }
        .set{ haplotype_table }

    haplotype_fasta
        .splitFasta(record: [id: true, seqString: true])
        .map{ [
            ['id': it[0].id, 'single_end': null, 'strandedness': null],
            it[1].id,
            it[1].seqString
        ] }
        .set{ haplotype_sequences }

    strain_table
        .combine(
            haplotype_table
                .join(haplotype_sequences, by: [0, 1])
                .mix(consensus_table),
            by: 0
        )
        .set{ answer_table }

    ANSWER_DUMP(answer_table)
}

process ANSWER_DUMP {
    input:
    tuple val(meta), file(bam), file(bai), val(accession), val(strain), val(name), val(freq), val(sequence)

    script:
    """
    echo "${meta.id}\t${accession}\t${strain}\t${name}\t${freq}" > ${name}.tsv
    echo ">${name}\n${sequence}" > ${name}.fasta
    """
}
