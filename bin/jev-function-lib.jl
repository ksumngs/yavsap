# JEV Functions Library
# Usage:
#    include("jev-function-lib.jl")

# Package imports
using XAM
using BioAlignments

"""
    baseatreferenceposition(record::BAM.Record, pos::Int)

Get the base at reference position `pos` present in the sequence of `record`.
"""
function baseatreferenceposition(record::BAM.Record, pos::Int)
    return BAM.sequence(record)[ref2seq(BAM.alignment(record), pos)[1]]
end # function

"""
    baseatreferenceposition(record::SAM.Record, pos::Int)

Get the base at reference position `pos` present in the sequence of `record`.
"""
function baseatreferenceposition(record::SAM.Record, pos::Int)
    return SAM.sequence(record)[ref2seq(SAM.alignment(record), pos)[1]]
end # function
