# JEV Functions Library
# Usage:
#    include("jev-function-lib.jl")

# Package imports
using BioAlignments
using BioSequences
using BioSymbols
using DataFrames
using Distributions
using FASTX
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
        return DNA_N
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

function Variant(vardict::Dict{String,Any})
    region  = vardict["chromosome"]
    pos     = vardict["position"]
    refbase = vardict["referencebase"]
    altbase = vardict["alternatebase"]
    tdepth  = vardict["totaldepth"]
    altdep  = vardict["alternatedepth"]

    refseq = LongDNASeq(refbase)
    altseq = LongDNASeq(altbase)

    return Variant(region, pos, refseq, altseq, tdepth, altdep)
end #function

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

function Haplotype(hapdict::Dict{String,Any})
    v = Variant.(hapdict["mutations"])
    return Haplotype(v)
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

function bamcounts2bamdata(bamcounts::AbstractVector{String})
    # Declare an empty bam stats data frame
    countsdata = DataFrame(
        chr                                  = String[],
        position                             = Int[],
        reference_base                       = String[],
        depth                                = Int[],
        base                                 = String[],
        count                                = Int[],
        avg_mapping_quality                  = Float64[],
        avg_basequality                      = Float64[],
        avg_se_mapping_quality               = Float64[],
        num_plus_strand                      = Int[],
        num_minus_strand                     = Int[],
        avg_pos_as_fraction                  = Float64[],
        avg_num_mismatches_as_fraction       = Float64[],
        avg_sum_mismatch_qualities           = Float64[],
        num_q2_containing_reads              = Int[],
        avg_distance_to_q2_start_in_q2_reads = Float64[],
        avg_clipped_length                   = Float64[],
        avg_distance_to_effective_3p_end     = Float64[]
    )

    # Transform the bam stats file
    for bamline in bamcounts
        # Split the base-independent stats by tabs
        bamfields = split(bamline, "\t")

        # Loop through the base-dependent stat blocks
        for i in 6:length(bamfields)
                # Split the base-dependent stats by colons
                basestats = split(bamfields[i], ":")

                # Parse the data into the correct types
                chr                                  = bamfields[1]
                position                             = parse(Int, bamfields[2])
                reference_base                       = bamfields[3]
                depth                                = parse(Int, bamfields[4])
                base                                 = basestats[1]
                count                                = parse(Int, basestats[2])
                avg_mapping_quality                  = parse(Float64, basestats[3])
                avg_basequality                      = parse(Float64, basestats[4])
                avg_se_mapping_quality               = parse(Float64, basestats[5])
                num_plus_strand                      = parse(Int, basestats[6])
                num_minus_strand                     = parse(Int, basestats[7])
                avg_pos_as_fraction                  = parse(Float64, basestats[8])
                avg_num_mismatches_as_fraction       = parse(Float64, basestats[9])
                avg_sum_mismatch_qualities           = parse(Float64, basestats[10])
                num_q2_containing_reads              = parse(Int, basestats[11])
                avg_distance_to_q2_start_in_q2_reads = parse(Float64, basestats[12])
                avg_clipped_length                   = parse(Float64, basestats[13])
                avg_distance_to_effective_3p_end     = parse(Float64, basestats[14])

                # Append the data to the dataframe
                push!(countsdata, [
                    chr,
                    position,
                    reference_base,
                    depth,
                    base,
                    count,
                    avg_mapping_quality,
                    avg_basequality,
                    avg_se_mapping_quality,
                    num_plus_strand,
                    num_minus_strand,
                    avg_pos_as_fraction,
                    avg_num_mismatches_as_fraction,
                    avg_sum_mismatch_qualities,
                    num_q2_containing_reads,
                    avg_distance_to_q2_start_in_q2_reads,
                    avg_clipped_length,
                    avg_distance_to_effective_3p_end
                ])
        end
    end

    return countsdata
end #function

function variantwarning(variant::DataFrameRow, reason::String)
    @warn string("Dropping variant '", variant.ALT, "' from ", variant.REGION, ":", variant.POS, "-", variant.POS, ". Reason: ", reason)
end #function

function matchingbam(variant::DataFrameRow, countsdata::DataFrame)
    # Find the matching bam stats for this variant call
    bammatches = countsdata[(countsdata.chr .== variant.REGION) .& (countsdata.position .== variant.POS) .& (countsdata.base .== variant.ALT),:]
    if length(eachrow(bammatches)) < 1
        variantwarning(variant, "variant mismatch")
        return nothing
    end
    bammatch = first(bammatches)

    # Sanity-check => ivar.REF == bam.reference_base
    if variant.REF != bammatch.reference_base
        error(string("iVar reports reference base ", variant.REF, ", but bam reports ", bammatch.reference_base, " at ", variant.REGION, ":", variant.POS))
    end

    return bammatch
end #function

function isdistributed(variant::DataFrameRow, countsdata::DataFrame; strandpos=0.1)
    bammatch = matchingbam(variant, countsdata)
    if isnothing(bammatch)
        return false
    end

    # Check that average read position is greater than 0.1 (variants are not called exclusively within 10% of the read edges)
    if abs(bammatch.avg_pos_as_fraction) < strandpos
        variantwarning(variant, string("read position ", bammatch.avg_pos_as_fraction))
        return false
    end

    return true
end #function

function isunbiased(variant::DataFrameRow, countsdata::DataFrame; maxdiff=0.1, maxskew=0.2, minp=1e-5)
    bammatch = matchingbam(variant, countsdata)
    if isnothing(bammatch)
        return false
    end

    # Pull the reference stats
    bamref = first(countsdata[(countsdata.chr .== variant.REGION) .& (countsdata.position .== variant.POS) .& (countsdata.base .== variant.REF),:])

    # Pull the contingencies into variables for clarity
    ref_plus  = bamref.num_plus_strand
    ref_minus = bamref.num_minus_strand
    alt_plus  = bammatch.num_plus_strand
    alt_minus = bammatch.num_minus_strand

    # If all reads are only on one strand (e.g. Nanopore), then skip this analysis
    if (ref_plus + alt_plus) == 0 || (ref_minus + alt_minus == 0)
        return true
    end

    # Get the proportions
    ref_total = ref_plus + ref_minus
    alt_total = alt_plus + alt_minus
    ref_plus_frac = ref_plus / ref_total
    alt_plus_frac = alt_plus / alt_total

    # Calculate difference
    diff = abs(ref_plus_frac - alt_plus_frac)

    # Calculate skew
    skew = abs(0.5 - alt_plus_frac)

    # Calculate probabilities
    p_val = prop_test([ref_plus, alt_plus], [ref_total, alt_total])

    # Check if this record is in range
    if (diff > maxdiff) && (skew > 2maxskew) && (p_val < minp)
        variantwarning(variant, string("strand distribution reference plus strand ", ref_plus, "/", ref_total, " alternate plus strand ", alt_plus, "/", alt_total, " p-value ", p_val))
        return false
    end

    return true

end

function yamlize(h::Haplotype; reason::Union{String,Nothing}=nothing)
    return string(
        "---\n",
        isnothing(reason) ? "" : string("reason: ", reason, "\n"),
        "mutations: \n",
        yamlize.(h.mutations)...
    )
end

function yamlize(v::Variant)
    return string(
        "  - chromosome: ",
        v.region,
        "\n",
        "    position: ",
        string(v.position),
        "\n",
        "    referencebase: ",
        string(v.referencebase),
        "\n",
        "    alternatebase: ",
        string(v.alternatebase),
        "\n",
        "    totaldepth: ",
        string(v.totaldepth),
        "\n",
        "    alternatedepth: ",
        string(v.alternatedepth),
        "\n",
    )
end

function mutaterecord(record::FASTA.Record, haplotype::Haplotype)
    mutationpos = position.(mutations(haplotype))
    newid = string(identifier(record),"_mutated_pos_",join(mutationpos, "_"))
    newseq = mutatesequence(sequence(record), haplotype)
    newdesc = ""
    if FASTA.hasdescription(record)
        newdesc = string(description(record), ", mutated at position(s) ", join(mutationpos, ", ", ", and "))
    end

    return FASTA.Record(
        newid,
        newdesc,
        newseq
    )

end

function mutatesequence(seq::NucleotideSeq, haplotype::Haplotype)
    newseq = seq

    for var in mutations(haplotype)
        newseq[var.position] = var.alternatebase[1]
    end

    return newseq
end
