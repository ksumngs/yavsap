include { HAPLOTYPE_YAML2TSV } from '../modules/local/haplotype-yaml2tsv.nf'
include { MINIMAP2_ALIGN } from '../modules/nf-core/modules/minimap2/align/main.nf'

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
    tree
    multiqc
    krona

    main:
    versions = Channel.empty()

    HAPLOTYPE_YAML2TSV(haplotype_yaml)
    HAPLOTYPE_YAML2TSV.out.tsv
        .map{ it[1] }
        .collectFile(name: 'haplotypes.tsv')
        .set{ haplotype_tsv }

    versions = versions.mix(HAPLOTYPE_YAML2TSV.out.versions)

    consensus_fasta
        .mix(haplotype_fasta)
        .map{ it[1] }
        .collectFile(name: 'haplotypes.fasta')
        .map{ [ [ 'id': 'haplotypes', 'single_end': true ], it ] }
        .set{ haplotype_sequences }

    MINIMAP2_ALIGN(haplotype_sequences, reference_fasta)

    MINIMAP2_ALIGN.out.paf.map{ it[1] }.set{ haplotype_alignment }

    versions = versions.mix(MINIMAP2_ALIGN.out.versions)

    accession
        .join(strain)
        .map{ [
            it[0].id,
            it[1],
            it[2]
        ] }
        .combine(haplotype_tsv.splitCsv(sep: '\t'), by: 0)
        .set{ strain_table }

    ECHO2TSV(strain_table)

    ECHO2TSV.out
        .collectFile(name: 'haplotype_strains.tsv')
        .set{ haplotype_strains }

    SEQUENCETABLE(
        haplotype_strains,
        haplotype_alignment,
        reference_fasta,
        tree,
        multiqc,
        krona
    )
}

process ECHO2TSV {
    label 'process_low'

    input:
    tuple val(sample), val(accession), val(strain), val(name), val(frequency)

    output:
    path "*.tsv"

    script:
    """
    echo "${sample}\t${accession}\t${strain}\t${name}\t${frequency}" \\
        > ${sample}_${accession}.tsv
    """
}

process YAVSAP_SUMMARY {
    input:
    file(sam)
    file(reference)
    file(tsv)

    script:
    """
    echo \$(ls)
    """
}
