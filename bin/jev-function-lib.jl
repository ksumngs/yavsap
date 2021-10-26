# JEV Functions Library
# Usage:
#    include("jev-function-lib.jl")

# Package imports
using BioAlignments
using BioSequences
using BioSymbols
using Combinatorics
using DataFrames
using Dates
using Distributions
using FASTX
using SHA
using StructArrays
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
    counts::AbstractArray{Int}
end

function haplotypecountsmutations(h::HaplotypeCounts)
    return h.haplotype.mutations
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

"""
    findoccurrences(haplotype::Haplotype, reads::Union{AbstractVector{BAM.Record}, BAM.Reader})

Finds the number of times variants occur within `reads` in accordance with `haplotype`.

Return an ``N``-dimensional array where the ``N``th dimension represents the ``N``th
variant call position within a haplotype. Position 1 of ``N`` indicates the number of times
the reference base was called for this variant position, while position 2 indicates the
number of times the alternate base was called. For example,
`findoccurences(hap, read)[1,2,1]` would give the number of times the reference was called
in position 1, the alternate was also called in position 2, and the reference was also
called in position 3.

"""
function findoccurrences(haplotype::Haplotype, reads::AbstractVector{BAM.Record})
    # Make things easier to call
    mutations = StructArray(haplotype.mutations)

    # Get only the reads that contain all of the variant positions
    containingreads = filter(b -> BAM.position(b) < min(mutations.position...) && BAM.rightposition(b) > max(mutations.position...), reads)

    # Set up haplotype counts
    hapcounts = zeros(Int, repeat([2], length(mutations))...)

    # Check every NGS read that contains both positions
    for record in containingreads
        # Pull the basecalls
        basecalls = baseatreferenceposition.([record], mutations.position)

        # Find how the basecalls stack up
        matches = matchvariant.(basecalls, Vector(mutations))

        if !any(matches .== :other)
            coordinate = CartesianIndex((Int.(matches .== :alternate) .+ 1)...)
            hapcounts[coordinate] += 1
        end #if
    end #for

    # Write the haplotype counts to the overall dataframe
    return HaplotypeCounts(haplotype, hapcounts)

end

function findoccurrences(haplotype::Haplotype, reads::BAM.Reader)
    return findoccurrences(haplotype, collect(reads))
end

function findsimulatedoccurrences(haplotype::Haplotype, reads::AbstractVector{BAM.Record}; iterations::Int=1000)
    # Extract the SNPs we care about
    mutations = haplotype.mutations

    # Create an empty array for the simulated long reads
    pseudoreads = Array{Symbol}(undef, iterations, length(mutations))

    # Set up haplotype counts
    hapcounts = zeros(Int, repeat([2], length(mutations))...)

    # Start iterating
    for i ∈ 1:iterations
        # Get the reads that contain the first mutation
        lastcontainingreads = filter(b -> BAM.position(b) < mutations[1].position && BAM.rightposition(b) > mutations[1].position, reads)

        # Pull a random read from that pool
        lastread = rand(lastcontainingreads)

        # Find this read's basecall at that position
        basecall = baseatreferenceposition(lastread, mutations[1].position)
        basematch = matchvariant(basecall, mutations[1])

        pseudoreads[i, 1] = basematch

        for j ∈ 2:length(mutations)
            if (BAM.position(lastread) < mutations[j].position && BAM.rightposition(lastread) > mutations[j].position)
                thisread = lastread
            else
                thiscontainingreads = filter(
                    b -> BAM.position(b) > BAM.rightposition(lastread) && BAM.position(b) < mutations[j].position && BAM.rightposition(b) > mutations[j].position,
                    reads
                )
                if length(thiscontainingreads) < 1
                    pseudoreads[i,j] = :other
                    continue
                end
                thisread = rand(thiscontainingreads)
            end #if

            # Find this read's basecall at that position
            basecall = baseatreferenceposition(thisread, mutations[j].position)
            basematch = matchvariant(basecall, mutations[j])

            pseudoreads[i, j] = basematch

            lastread = thisread

        end #for

        matches = pseudoreads[i, :]
        if !any(matches .== :other)
            coordinate = CartesianIndex((Int.(matches .== :alternate) .+ 1)...)
            hapcounts[coordinate] += 1
        end #if

    end #for

    return HaplotypeCounts(haplotype, hapcounts)

end #function

"""
    linkage(counts::AbstractArray{Int})

Calculates the linkage disequilibrium and Chi-squared significance level of a combination of
haplotypes whose number of occurrences are given by `counts`.

`counts` is an ``N``-dimensional array where the ``N``th dimension represents the ``N``th
variant call position within a haplotype. `findoccurrences` produces such an array.
"""
function linkage(counts::AbstractArray{Int})
    # Get the probability of finding a perfect reference sequence
    P_allref = first(counts) / sum(counts)

    # Get the probabilities of finding reference bases in any of the haplotypes
    P_refs = sumsliced.([counts], 1:ndims(counts)) ./ sum(counts)

    # Calculate linkage disequilibrium
    Δ = P_allref - prod(P_refs)

    # Calculate the test statistic
    r = Δ / (prod(P_refs .* (1 .- P_refs))^(1/ndims(counts)))
    Χ_squared = r^2 * sum(counts)

    # Calculate the significance
    p = 1 - cdf(Chisq(1), Χ_squared)

    return Δ, p
end #function

"""
    sumsliced(A::AbstractArray, dim::Int, pos::Int=1)

Sum all elements that are that can be referenced by `pos` in the `dim` dimension of `A`.

# Example

```julia-repl
julia> A = reshape(1:8, 2, 2, 2)
2×2×2 reshape(::UnitRange{Int64}, 2, 2, 2) with eltype Int64:
[:, :, 1] =
 1  3
 2  4

[:, :, 2] =
 5  7
 6  8

julia> sumsliced(A, 2)
16

julia> sumsliced(A, 2, 2)
20
```

Heavily inspired by Holy, Tim "Multidimensional algorithms and iteration"
<https://julialang.org/blog/2016/02/iteration/#filtering_along_a_specified_dimension_exploiting_multiple_indexes>
"""
function sumsliced(A::AbstractArray, dim::Int, pos::Int=1)
    i_pre  = CartesianIndices(size(A)[1:dim-1])
    i_post = CartesianIndices(size(A)[dim+1:end])
    return sum(A[i_pre, pos, i_post])
end

function writesimulatedhaplotypes(variants, bamreads, outyaml, max_var, hap_min, hap_sig)
    # Hold any variants that have already been shown to exhibit linkage
    linkedvariants = Variant[]

    @info string("# of variants:", '\t', length(variants))

    # Find the number of variants allowed to be in a haplotype
    max_var = min(max_var, length(variants))

    # Find all possible pairings of these variants
    for i ∈ max_var:-1:2
        @info string("Now considering haplotypes of ", '\t', i, '\t', "variants.", '\t', Dates.now())

        # Take out any variants that have already been considered
        filteredvariants = filter(!v -> v ∈ linkedvariants, variants)

        # Find the combinations of variants that can be assembled into haplotypes
        variantcombos = combinations(filteredvariants, i)
        @info string('\t', length(variantcombos), '\t', "haplotype(s) possible")

        # Set up a rough progressbar
        print(string('\t', '\t', "progress: "))
        numcombos = length(variantcombos)
        j = 0
        lastpoint = 0
        endpoint = 57

        # Compute for each possible combo
        for variantcombo in variantcombos
            j += 1
            percentcomplete = j / numcombos
            if percentcomplete > lastpoint / endpoint
                lastpoint += 1
                print(".")
            end
            # Convert the combo into a haplotype
            expectedhaplotype = Haplotype(variantcombo)

            # Find the simulated counts
            haplotypecounts = findsimulatedoccurrences(expectedhaplotype, bamreads)

            # Check if this haplotype is significant
            linkstats = linkage(haplotypecounts.counts)
            if linkstats[2] <= hap_sig && last(haplotypecounts.counts) >= hap_min

                # Add the variant positions of this haplotype to the list of variants that can
                # be ignored
                push!(linkedvariants, cat(haplotypecounts.haplotype.mutations, dims=1)...)

                # Write this haplotype to the output file
                open(outyaml, "a") do f
                    write(f, string(yamlize(haplotypecounts, reason="linkage")))
                end #do
            end #if
        end #for
        print("\n")
    end #for

    return linkedvariants
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

function yamlize(hc::HaplotypeCounts; reason::Union{String,Nothing}=nothing)
    occurrences = "occurrences:\n"
    for i in CartesianIndices(hc.counts)
        location = [Tuple(i)...]
        variantpattern = string.(replace(replace(location, 1 => "ref"), 2 => "alt"))
        key = join(variantpattern, "_")
        occurrences = string(occurrences, "  ", key, ": ", hc.counts[i], "\n")
    end

    return string(
        yamlize(hc.haplotype, reason=reason),
        occurrences,
        "ldcofficient: ",
        linkage(hc.counts)[1],
        "\n",
        "p_val: ",
        linkage(hc.counts)[2],
        "\n"
    )
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
    newseq = mutatesequence(sequence(record), haplotype)
    sequencehash = bytes2hex(sha1(string(newseq)))
    newid = sequencehash[1:8]
    newdesc = string(description(record), ", variant ", sequencehash)

    return FASTA.Record(
        newid,
        newdesc,
        newseq
    )

end

function mutatesequence(seq::NucleotideSeq, haplotype::Haplotype)
    newseq = seq

    for var in haplotype.mutations
        newseq[var.position] = var.alternatebase[1]
    end

    return newseq
end
