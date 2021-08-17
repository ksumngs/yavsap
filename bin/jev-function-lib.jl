# JEV Functions Library
# Usage:
#    include("jev-function-lib.jl")

# Package imports
using BioAlignments
using BioSequences
using BioSymbols
using DataFrames
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
