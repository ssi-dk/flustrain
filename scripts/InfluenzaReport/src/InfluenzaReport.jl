"""
    InfluenzaReport

This package is the Julia-side (backend) of the flupipe Snakemake pipeline.
It's purpose is to be *application-specific* to the flupipe, NOT to include
generally useful functions. For those, see `Influenza.`
"""
module InfluenzaReport

using FASTX
using BioSequences
using BioAlignments
using CodecZlib
using ErrorTypes
using Influenza
using Printf
using Serialization
using Transducers
using Folds
using Plots

const N_SEGMENTS = length(instances(Segment))
const SegmentTuple{T} = NTuple{N_SEGMENTS, T}
const TERMINAL = 25

function split!(v::Vector{SubString{String}}, s::Union{String, SubString{String}}, sep::UInt8)
    n = 0
    start = 1
    @inbounds for i in 1:ncodeunits(s)
        if codeunit(s, i) == sep
            n += 1
            n >= length(v) && throw(BoundsError(v, n+1))
            substr = SubString(s, start, i-1)
            v[n] = substr
            start = i + 1
        end
    end
    @inbounds v[n+1] = SubString(s, start, ncodeunits(s))
    v
end

@enum Severity::UInt8 trivial important

struct ErrorMessage
    severity::Severity
    msg::String
end

is_trivial(x::ErrorMessage) = x.severity == trivial

function format(msg::ErrorMessage)
    preface = msg.severity == trivial ? "         " :
                                        "ERROR:   "
    preface * msg.msg
end

const _IMPORTANT = Tuple(Bool[1,1,0,0,1,0,1,1,1,1,1,1,1,0,1,1,1])
@assert length(_IMPORTANT) == length(instances(Protein))
is_important(x::Protein) = @inbounds _IMPORTANT[reinterpret(UInt8, x) + 0x01]

include("assembly.jl")
include("depths.jl")
include("reference.jl")

include("alignment.jl")
include("report.jl")


# Actually call the stuff here (else, we import "interactive.jl")
if abspath(PROGRAM_FILE) == @__FILE__
    # TODO: Add GridIon functionality
    if length(ARGS) != 2
        println("Usage: julia InfluenzaReport.jl out_dir ref_dir")
        exit(1)
    end

    outdir = first(ARGS)
    refdir = last(ARGS)

    reportpath = joinpath(outdir, "report.txt")
    alndir = joinpath(outdir, "aln")
    consdir = joinpath(outdir, "consensus")
    plotdir = joinpath(outdir, "depths")

    illumina_snakemake_entrypoint(reportpath, refdir, alndir, consdir, plotdir)
else
    include("interactive.jl")
end

end # module
