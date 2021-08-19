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

struct Haplotype
    mutations::AbstractVector{Variant}
end

function Haplotype(var::Variant)
    return Haplotype([var])
end

function Base.position(var::Variant)
    return var.position
end

struct HaplotypeCounts
    haplotype::Haplotype
    counts::AbstractMatrix{Int}
end

function HaplotypeCounts(haplotype::Haplotype, counts::Int...)
    if !isinteger(sqrt(length(counts)))
        error(ArgumentError("number of count arguments must be reshabable to a square matrix"))
    end
    nhaps = Int(sqrt(length(counts)))
    if length(haplotype.mutations) != nhaps
        error(ArgumentError("the number of count arguments should be equal to the square of the number of variants in the haplotype"))
    end
    countmat = Matrix(reshape(counts, (nhaps, nhaps)))
    return HaplotypeCounts(haplotype, countmat)
end

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

"""
    matchvariant(base::Union{NucleotideSeq,DNA,AbstractVector{DNA}}, var::Variant)

Checks if `base` matches the reference or variant expected in `var`, and returns a symbol
indicating which, if any, it matches.

Returned values can be `:reference` for a reference match, `:alternate` for an alternate
match, or `:other` for no match with the given variant.
"""
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

function findoccurrences(haplotype::Haplotype, reads::AbstractVector{BAM.Record})
    # Set up haplotype counts
    ref_ref = 0
    ref_alt = 0
    alt_ref = 0
    alt_alt = 0

    # Make things easier to call
    mutations = haplotype.mutations

    # Get only the reads that contain all of the variant positions
    containingreads = filter(b -> BAM.position(b) < min(position.(mutations)...) && BAM.rightposition(b) > max(position.(mutations)...), reads)

    # Check every NGS read that contains both positions
    for record in containingreads
        # Find this read's basecalls for both locations
        try
            basecalls = baseatreferenceposition.([record], position.(mutations))

            # Find how the basecalls stack up
            matches = matchvariant.(basecalls, mutations)

            # Find which haplotype we're dealing with, if any
            # Theoretically, we could make an n-dimensional array for haplotypes consisting
            # of n mutations, but for now we'll stick with two
            if     matches[1] == :reference && matches[2] == :reference
                ref_ref += 1
            elseif matches[1] == :reference && matches[2] == :alternate
                ref_alt += 1
            elseif matches[1] == :alternate && matches[2] == :reference
                alt_ref += 1
            elseif matches[1] == :alternate && matches[2] == :alternate
                alt_alt += 1
            end #if
        catch e
            @warn "Base out of bounds" record e
        end #try
    end #for

    # Write the haplotype counts to the overall dataframe
    return HaplotypeCounts(haplotype, [ref_ref ref_alt; alt_ref alt_alt])

end

"""
    haplotypeoccurances(variantparings::Vector{Vector{Variant}}, reader::BAM.Reader)

Tabulates the number of types each haplotype corresponding to the reference/alternate
pairings present in `variantpairings` occurs in the reads contained in `reader`.
"""
function haplotypeoccurances(variantparings::Vector{Vector{Variant}}, reader::BAM.Reader)
    # Declare an empty dataframe for the possible haplotypes
    haplotypedata = DataFrame(
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
    for varpair in variantparings
        # Set up haplotype counts
        ref_ref = 0
        ref_alt = 0
        alt_ref = 0
        alt_alt = 0

        # Get the reads that contain the first variant position
        containingreads = collect(eachoverlap(reader, varpair[1].region, varpair[1].position:varpair[1].position))

        # Get the reads that contain the second variant position
        containingreads = containingreads[BAM.rightposition.(containingreads) .>= varpair[2].position]

        # Check every NGS read that contains both positions
        for record in containingreads
            # Find this read's basecalls for both locations
            try
                base1 = baseatreferenceposition(record, varpair[1].position)
                base2 = baseatreferenceposition(record, varpair[2].position)

                # Find how the basecalls stack up
                match1 = matchvariant(base1, varpair[1])
                match2 = matchvariant(base2, varpair[2])

                # Find which haplotype we're dealing with, if any
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

        # Write the haplotype counts to the overall dataframe
        if sum([ref_ref, alt_ref, ref_alt, alt_alt]) > 0
            push!(haplotypedata, [varpair[1].position, varpair[2].position, varpair[1].alternatebase, varpair[2].alternatebase, ref_ref, ref_alt, alt_ref, alt_alt])
        end #if
    end #for

    return haplotypedata
end #function

"""
    appendlinkagestatistics!(haplotypecounts::DataFrame)

Calculates linkage disequilibrium and the Χ-squared p-value of linkage disequilibrium for
haplotypes contained in `haplotypecounts`, filtering out rows where it is invalid to
calculate linkage disequilibrium. `haplotypecounts` is expected to be in the same format
returned by [`haplotypeoccurances`](@ref).
"""
function appendlinkagestatistics!(haplotypecounts::DataFrame)
    # Filter out problem cases
    # 1. Filter out instances where linkage disequilibrium can't be calculated (0 instances of an allele combo)
    filter!(v -> prod([v.Reference_Reference, v.Reference_Alternate, v.Alternate_Reference, v.Alternate_Alternate]) > 0, haplotypecounts)
    # 2. Filter out low read depth instances (< 100 total reads)
    filter!(v -> sum([v.Reference_Reference, v.Reference_Alternate, v.Alternate_Reference, v.Alternate_Alternate]) < 100, haplotypecounts)
    # 3. Filter out low double-mutation instances (< 5 alt_alt reads)
    filter!(v -> v.Alternate_Alternate < 5, haplotypecounts)

    # Split out the fields into managable variables
    ref_ref = haplotypecounts.Reference_Reference
    ref_alt = haplotypecounts.Reference_Alternate
    alt_ref = haplotypecounts.Alternate_Reference
    alt_alt = haplotypecounts.Alternate_Alternate

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
    haplotypecounts.Linkage_Disequilibrium = linkage
    haplotypecounts.Linkage_Significance = sig

    return haplotypecounts
end
