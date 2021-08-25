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
    seqpos = myref2seq(BAM.alignment(record), pos)[1]
    if seqpos > 0 && seqpos < BAM.seqlength(record)
        return BAM.sequence(record)[seqpos]
    else
        return dna"N"
    end
end # function

struct Variant
    region::String
    position::Int
    referencebase::NucleotideSeq
    alternatebase::NucleotideSeq
    totaldepth::Int
    alternatedepth::Int
end

function Variant(data::DataFrameRow)
    region = data.REGION
    pos = data.POS
    refbase = data.REF
    altbase = data.ALT
    tdepth = data.TOTAL_DP
    altdep = data.ALT_DP

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

    return Variant(region, pos, refseq, altseq, tdepth, altdep)

end # function

function region(v::Variant)
    return v.region
end

function Base.position(v::Variant)
    return v.position
end

function referencebase(v::Variant)
    return v.referencebase
end

function alternatebase(v::Variant)
    return v.alternatebase
end

function frequency(v::Variant)
    return v.alternatedepth / v.totaldepth
end

function totaldepth(v::Variant)
    return v.totaldepth
end

function alternatedepth(v::Variant)
    return v.alternatedepth
end

struct Haplotype
    mutations::AbstractVector{Variant}
end

function Haplotype(var::Variant)
    return Haplotype([var])
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

function counts(h::HaplotypeCounts)
    return h.counts
end

function haplotype(h::HaplotypeCounts)
    return h.haplotype
end

function mutations(h::HaplotypeCounts)
    return h.haplotype.mutations
end

function mutations(h::Haplotype)
    return h.mutations
end

function ref_ref(h::HaplotypeCounts)
    return h.counts[1,1]
end

function ref_alt(h::HaplotypeCounts)
    return h.counts[1,2]
end

function alt_ref(h::HaplotypeCounts)
    return h.counts[2,1]
end

function alt_alt(h::HaplotypeCounts)
    return h.counts[2,2]
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

    # Check that the requested base is in range
    if !seqisinrange(workingalignment, i)
        return (0, OP_HARD_CLIP)
    end

    # Perform regular alignment search, minus any hard clipping
    return ref2seq(workingalignment, i)

end #function

function seqisinrange(aln::Alignment, i::Int)
    reflen = i - first(aln.anchors).refpos
    seqlen = last(aln.anchors).seqpos - first(aln.anchors).seqpos
    return seqlen > reflen
end #function

function firstseqpos(aln::Alignment)
    return first(aln.anchors).seqpos
end #function

function lastseqpos(aln::Alignment)
    return last(aln.anchors).seqpos
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
    # Make things easier to call
    mutations = haplotype.mutations

    # Get only the reads that contain all of the variant positions
    containingreads = filter(b -> BAM.position(b) < min(position.(mutations)...) && BAM.rightposition(b) > max(position.(mutations)...), reads)

    # Set up haplotype counts
    hapcounts = zeros(Int, (length(containingreads), 4))

    # Check every NGS read that contains both positions
    Threads.@threads for i in 1:length(containingreads)
        # Extract the pertinant record (like this was a for each loop)
        record = containingreads[i]

        # Pull the basecalls
        basecalls = baseatreferenceposition.([record], position.(mutations))

        # Find how the basecalls stack up
        matches = matchvariant.(basecalls, mutations)

        # Find which haplotype we're dealing with, if any
        # Theoretically, we could make an n-dimensional array for haplotypes consisting
        # of n mutations, but for now we'll stick with two
        if     matches[1] == :reference && matches[2] == :reference
            hapcounts[i,1] = 1
        elseif matches[1] == :reference && matches[2] == :alternate
            hapcounts[i,2] = 1
        elseif matches[1] == :alternate && matches[2] == :reference
            hapcounts[i,3] = 1
        elseif matches[1] == :alternate && matches[2] == :alternate
            hapcounts[i,4] = 1
        end #if
    end #for

    # Tabulate results
    ref_ref = sum(hapcounts[:,1])
    ref_alt = sum(hapcounts[:,2])
    alt_ref = sum(hapcounts[:,3])
    alt_alt = sum(hapcounts[:,4])

    # Write the haplotype counts to the overall dataframe
    return HaplotypeCounts(haplotype, [ref_ref ref_alt; alt_ref alt_alt])

end

function findoccurrences(haplotype::Haplotype, reads::BAM.Reader)
    return findoccurrences(haplotype, collect(reads))
end

function linkage(counts::AbstractMatrix{Int})
    # Split out the fields into managable variables
    ref_ref = counts[1,1]
    ref_alt = counts[1,2]
    alt_ref = counts[2,1]
    alt_alt = counts[2,2]

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
    Δ = p_ref_ref .- (p_ref_1 .* p_ref_2)

    # Get the test statstic for linkage disequilibrium
    r = Δ / sqrt(p_ref_1 * (1 - p_ref_1) * p_ref_2 * (1 - p_ref_2))
    Χ_squared = r^2 * total

    # Check for significance
    p = 1 .- cdf.(Chisq(1), Χ_squared)

    return Δ, p
end #function

function prop_test(x, n)
    k = length(x)
    conf = 0.95
    estimate = x ./ n
    yates = 0.5
    Δ = estimate[1] - estimate[2]
    yates = min(yates, abs(Δ) / sum(1 ./ n))
    width = quantile(Normal(0,1), (1 + conf) / 2) * sqrt(sum(estimate .* (1 .- estimate) ./ n)) + yates * sum(1 ./ n)
    cint = [max(Δ-width, -1), min(Δ+width, 1)]
    p = sum(x) / sum(n)
    param = k - 1
    x2 = cat(x, n .- x, dims=2)
    E = cat(n .* p, n .* (1-p), dims=2)
    statistic = sum((abs.(x2 .- E) .- yates).^2 ./ E)
    p_val = 1 .- cdf.(Chisq(param), statistic)
end #function
