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

const SegmentTuple{T} = NTuple{length(instances(Segment)), T}
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

end # module
