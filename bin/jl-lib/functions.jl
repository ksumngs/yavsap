function ncbi_link(accession_number::AbstractString)
    return a(accession_number; href="https://ncbi.nlm.nih.gov/nuccore/$accession_number")
end #function

function base_table(sam::SAM.Record, refseq::NucleotideSeq)
    sequence = SAM.sequence(sam)
    alignment = SAM.alignment(sam)
    seq_start = SAM.position(sam)
    seq_end = SAM.rightposition(sam)

    basecells = []

    for (i, refbase) in enumerate(refseq)
        if i < seq_start || i > seq_end
            basetext = "-"
            baseclass = ""
        else
            (loc, op) = ref2seq(alignment, i)
            if isdeleteop(op) || ismetaop(op)
                basetext = "-"
                baseclass = "variant"
            else
                altbase = sequence[loc]
                basetext = string(altbase)
                baseclass = string(altbase)
                if altbase != refbase
                    baseclass = baseclass * " variant"
                end #if
            end #if
        end #if

        push!(basecells, td(basetext; class=baseclass))
    end #for

    return basecells
end #function

function sample_rows(
    samplename::AbstractString,
    reference::NucleotideSeq,
    alignments::AbstractVector{SAM.Record},
    data::AbstractDataFrame,
)
    sample_table = filter(v -> v.sample == samplename, data)
    sample_strain = replace(first(sample_table.strain), "_" => ": ")
    sample_accession = first(sample_table.accession)

    num_strains = length(eachrow(sample_table))

    # Get the alignment of this sample's consensus sequence
    consensus_record = first(
        filter(s -> contains(SAM.tempname(s), "Consensus_$(samplename)"), alignments)
    )

    rows = EzXML.Node[]
    push!(
        rows,
        tr(
            th(samplename; rowspan=(num_strains + 1)),
            td("Consensus"; colspan=2),
            td(sample_strain),
            td(ncbi_link(sample_accession)),
            base_table(consensus_record, reference)...,
        ),
    )

    for haplotype_table in eachrow(sample_table)
        haplotype_name = haplotype_table.haplotype
        haplotype_frequency = haplotype_table.frequency

        haplotype_record = first(
            filter(
                s -> contains(SAM.tempname(s), haplotype_table.haplotype), alignment_records
            ),
        )

        push!(
            rows,
            tr(
                td(haplotype_name),
                td(em("$(round(haplotype_frequency * 100))%")),
                td(sample_strain),
                td(ncbi_link(sample_accession)),
                base_table(haplotype_record, reference)...,
            ),
        )
    end #for

    return rows
end
