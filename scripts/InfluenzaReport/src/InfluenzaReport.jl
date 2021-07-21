"""
    InfluenzaReport

This package is the Julia-side (backend) of the flupipe Snakemake pipeline.
It's purpose is to be *application-specific* to the flupipe, NOT to include
generally useful functions. For those, see `Influenza.jl`
"""
module InfluenzaReport

using FASTX
using BioSequences
using ErrorTypes
using Influenza
using KMATools
using Printf
using Transducers
using Folds
using Plots
using CodecZlib

using Influenza: Assembly, Reference, AlignedAssembly

const N_SEGMENTS = length(instances(Segment))
const SegmentTuple{T} = NTuple{N_SEGMENTS, T}
const TERMINAL = 25

"Too many indel errors in sequence - alignment probably went wrong"
struct ErrorTooManyIndels <: Influenza.ProteinError
    n::UInt32
end

function Base.print(io::IO, x::ErrorTooManyIndels)
    print(io, "Too many indel errors, found ", x.n)
end

# TODO: Difference between depths in kma and kma2? Think about it!
# TODO: Superinfection?

const _IMPORTANT = Tuple(Bool[1,1,0,0,1,0,1,1,1,1,1,1,1,0,1,1,1])
@assert length(_IMPORTANT) == length(instances(Protein))
is_important(x::Protein) = @inbounds _IMPORTANT[reinterpret(UInt8, x) + 0x01]

include("alignedassembly.jl")
include("depths.jl")
include("report.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 3
        println("Usage: julia InfluenzaReport.jl platform out_dir ref_dir")
        exit(1)
    end
    platform, outdir, refdir = ARGS
    illumina = if platform == "illumina"
        true
    elseif platform == "nanopore"
        false
    else
        println("Platform must be \"illumina\" or \"nanopore\"")
        exit(1)
    end

    reportpath = joinpath(outdir, "report.txt")
    alndir = joinpath(outdir, "tmp/aln")
    consdir = joinpath(outdir, "consensus")
    plotdir = joinpath(outdir, "depths")
    
    if illumina
        illumina_snakemake_entrypoint(reportpath, refdir, alndir, consdir, plotdir)
    else
        nanopore_snakemake_entrypoint(reportpath, refdir, alndir, consdir, plotdir)
    end
end

end # module
