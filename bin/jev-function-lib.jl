# JEV Functions Library
# Usage:
#    include("jev-function-lib.jl")

# Package imports
using BioAlignments
using BioSequences
using BioSymbols
using DataFrames
using Distributions
using XAM

"""
    baseatreferenceposition(record::BAM.Record, pos::Int)

Get the base at reference position `pos` present in the sequence of `record`.
"""
function baseatreferenceposition(record::BAM.Record, pos::Int)
    return BAM.sequence(record)[myref2seq(BAM.alignment(record), pos)[1]]
end # function

"""
    baseatreferenceposition(record::SAM.Record, pos::Int)

Get the base at reference position `pos` present in the sequence of `record`.
"""
function baseatreferenceposition(record::SAM.Record, pos::Int)
    return SAM.sequence(record)[myref2seq(SAM.alignment(record), pos)[1]]
end # function

struct Variant
    region::String
    position::Int
    referencebase::NucleotideSeq
    alternatebase::NucleotideSeq
end

function Variant(data::DataFrameRow)
    region = data.REGION
    pos = data.POS
    refbase = data.REF
    altbase = data.ALT

    # Check for insertion
    if first(altbase) == '+'
        altbase = string(refbase, altbase[2:end])
    end # if

    # Check for deletion
    if first(altbase) == '-'
        altbase = "-"
    end

    refseq = LongDNASeq(refbase)
    altseq = LongDNASeq(altbase)

    return Variant(region, pos, refseq, altseq)

end # function

"""
    myref2seq(aln::Alignment, i::Int)

Replicates the functionality of BioAlignments `ref2seq`, but can handle hard clips
by effectively removing them for the intent of finding the position.
"""
function myref2seq(aln::Alignment, i::Int)
    if aln.anchors[2].op == OP_HARD_CLIP
        # Hard clipping was shown on operation 2
        # (operation 1 is always a start position)

        # Save where the clipping ends
        alnstart = aln.anchors[2]

        # Declare a new empty array where we can rebuild the alignment
        newanchors = AlignmentAnchor[]

        # Rebase the start of our new alignment to where the clipping ends
        push!(newanchors, AlignmentAnchor(
            0,
            aln.anchors[1].refpos,
            OP_START
        ))

        # Add new anchors
        for j in 3:(length(aln.anchors)-1)
            newanchor = AlignmentAnchor(
                aln.anchors[j].seqpos - alnstart.seqpos,
                aln.anchors[j].refpos,
                aln.anchors[j].op
            )
            push!(newanchors, newanchor)
        end #for

        # Package up our new alignment
        workingalignment = Alignment(newanchors)
    else
        # Package up the old alignment if there was no hard clipping
        workingalignment = aln
    end #if

    # Perform regular alignment search, minus any hard clipping
    return ref2seq(workingalignment, i)

end #function

function matchvariant(base::NucleotideSeq, var::Variant)
    refbase = LongDNASeq(var.referencebase)
    altbase = LongDNASeq(var.alternatebase)

    if base == refbase
        return :reference
    elseif base == altbase
        return :alternate
    else
        return :other
    end #if
end #function

function matchvariant(base::DNA, var::Variant)
    return matchvariant(LongDNASeq([base]), var)
end

function matchvariant(base::AbstractVector{DNA}, var::Variant)
    return matchvariant(LongDNASeq(base), var)
end

function findvariantlinkages(variantcombos::Vector{Vector{Variant}}, reader::BAM.Reader)
    # Declare an empty dataframe for the possible combinations
    variantcombodata = DataFrame(
        Position_1 = Int[],
        Position_2 = Int[],
        Call_1 = String[],
        Call_2 = String[],
        Reference_Reference = Int[],
        Reference_Alternate = Int[],
        Alternate_Reference = Int[],
        Alternate_Alternate = Int[]
    )


    # Find the reads that match every possible combo
    for variantpair in variantcombos
        # Set up variant combo counts
        ref_ref = 0
        ref_alt = 0
        alt_ref = 0
        alt_alt = 0

        # Get the reads that contain the first variant position
        containingreads = collect(eachoverlap(reader, variantpair[1].region, variantpair[1].position:variantpair[1].position))

        # Get the reads that contain the second variant position
        containingreads = containingreads[BAM.rightposition.(containingreads) .>= variantpair[2].position]

        # Check every NGS read that contains both positions
        for record in containingreads
            # Find this read's basecalls for both locations
            try
                base1 = baseatreferenceposition(record, variantpair[1].position)
                base2 = baseatreferenceposition(record, variantpair[2].position)

                # Find how the basecalls stack up
                match1 = matchvariant(base1, variantpair[1])
                match2 = matchvariant(base2, variantpair[2])

                # Find which combination we're dealing with, if any
                if     match1 == :reference && match2 == :reference
                    ref_ref += 1
                elseif match1 == :reference && match2 == :alternate
                    ref_alt += 1
                elseif match1 == :alternate && match2 == :reference
                    alt_ref += 1
                elseif match1 == :alternate && match2 == :alternate
                    alt_alt += 1
                end #if
            catch e
                @warn "Base out of bounds" record e
            end #try
        end #for

        # Write which combos we have to the overall dataframe
        if sum([ref_ref, alt_ref, ref_alt, alt_alt]) > 0
            push!(variantcombodata, [variantpair[1].position, variantpair[2].position, variantpair[1].alternatebase, variantpair[2].alternatebase, ref_ref, ref_alt, alt_ref, alt_alt])
        end #if
    end #for

    return variantcombodata
end #function

function appendlinkagestatistics!(var_combos::DataFrame)
    filter!(v -> prod([v.Reference_Reference, v.Reference_Alternate, v.Alternate_Reference, v.Alternate_Alternate]) > 0, var_combos)
    # Split out the fields into managable variables
    ref_ref = var_combos.Reference_Reference
    ref_alt = var_combos.Reference_Alternate
    alt_ref = var_combos.Alternate_Reference
    alt_alt = var_combos.Alternate_Alternate

    # Tabulate allele frequencies
    total = ref_ref .+ ref_alt .+ alt_ref .+ alt_alt
    ref_1 = ref_ref .+ ref_alt
    ref_2 = ref_ref .+ alt_ref
    alt_1 = alt_ref .+ alt_alt
    alt_2 = ref_alt .+ alt_alt

    # Calculate allele probabilities (unused figures are not present)
    p_ref_ref = ref_ref ./ total
    p_ref_1 = ref_1 ./ total
    p_ref_2 = ref_2 ./ total

    # Calculate linkage disequilibrium
    linkage =
        p_ref_ref .- (p_ref_1 .* p_ref_2)

    # Get the test statstic for linkage disequilibrium
    Δlinkage = (sqrt.(total .* alt_alt) .- sqrt.(alt_1 .* alt_2)) ./ total
    linkageerror = (0.5 ./ total) .* sqrt.(
        ref_ref .- (4 .* total .* Δlinkage .*
            (
                ((alt_1 .+ alt_2)./(2 .* sqrt.(alt_1 .* alt_2))) .-
                (sqrt.(alt_1 .* alt_2) ./ total)
            )
        )
    )
    linkagestatistic = abs.(Δlinkage ./ linkageerror)

    # Check for significance
    sig = 1 .- cdf.(Chisq(1), linkagestatistic)

    # Add linkage and significance to the dataframe
    var_combos.Linkage_Disequilibrium = linkage
    var_combos.Linkage_Significance = sig

    return var_combos
end
