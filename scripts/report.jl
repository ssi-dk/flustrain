using FASTX
using BioSequences
using CodecZlib

const TERMINAL = 25

function get_depth(path::AbstractString)
    reference = DNA[]
    depths = UInt32[]

    open(path) do io
        lines = eachline(GzipDecompressorStream(io))
        iterate(lines) # skip header
        stripped = (rstrip(line) for line in lines)
        for fields in (split(strip) for strip in stripped if !isempty(strip))
            nucleotide = convert(DNA, first(first(fields)))
            nucleotide == DNA_Gap && continue

            # A, C, G, T and also -, but not N
            depth = sum(x -> parse(UInt32, x), @view fields[2:5]) + parse(UInt32, fields[7])
            push!(depths, depth)
            push!(reference, nucleotide)
        end
        (depths, LongDNASeq(reference))
    end
end

# Depth naturally declines towards the ends of the template. This is not
# necessarily a sign of bad sequencing. So I include the first and last 25 bp,
# and check for depth < 25.
count_low_depths(x::Vector{<:Integer}) = count(<(25), @view x[1+TERMINAL:end-TERMINAL])

function count_insignificant_bases(rec::FASTA.Record)
    count(byte -> islowercase(Char(byte)), @view rec.data[rec.sequence])
end

count_ns(rec::FASTA.Record) = count(@view rec.data[rec.sequence]) do byte
    (byte === UInt8('n')) | (byte === UInt8('N'))
end

struct ConsensusState
    segment::String
    ok::Bool
    mean_depth::Float64
    coverage::Float64
    n_low_depth::Int
    n_insignificant::Int
    n_ns::Int
end

function ConsensusState(segment::AbstractString, conspath::AbstractString,
    matpath::AbstractString
)::ConsensusState
    # If file is empty
    first_it = open(iterate, FASTA.Reader, conspath)
    EMPTY_STATE = ConsensusState(segment, false, 0.0, 0.0, 0, 0, 0)
    isnothing(first_it) && return EMPTY_STATE

    consensus::FASTA.Record = first(first_it)
    isempty(consensus.sequence) && return EMPTY_STATE

    depths, reference = get_depth(matpath)
    n_low_depth = count_low_depths(depths)
    n_ns = count_ns(consensus)
    n_insignificant = count_insignificant_bases(consensus)
    mean_depth = sum(depths) / length(depths)
    coverage = count(!iszero, depths) / length(depths)

    ok = iszero(n_low_depth) && iszero(n_ns) && iszero(n_insignificant)

    return ConsensusState(segment, ok, mean_depth, coverage, n_low_depth, n_insignificant, n_ns)
end

function report_text(state::ConsensusState)::String
    dep = round(state.mean_depth, digits=2)
    cov = round(state.coverage, digits=2)
    lines = ["$(state.segment):\tDepth: $dep\tCoverage: $cov"]
    # If the state is OK, we just return that.
    state.ok && return first(lines)

    if !iszero(state.n_low_depth)
        push!(lines, "\t$(state.n_low_depth) bases (excluding terminal $(TERMINAL)bp) has depth < 10")
    end
    if !iszero(state.n_ns)
        push!(lines, "\t$(state.n_ns) bases called as \"N\"")
    end
    if !iszero(state.n_insignificant)
        push!(lines, "\t$(state.n_insignificant) basecalls without statistical significance")
    end
    join(lines, '\n')
end

function main(conspaths::Vector{<:AbstractString}, matpaths::Vector{<:AbstractString},
    accpath::AbstractString, reportpath::AbstractString
)
    # Step one: Calculate the report info for each segment
    states = ConsensusState[]
    for (conspath, matpath) in zip(conspaths, matpaths)
        segment = first(splitext(basename(conspath)))
        push!(states, ConsensusState(segment, conspath, matpath))
    end

    # Then print the accepted ones:
    open(accpath, "w") do accfile
        for state in states
            state.ok && println(accfile, state.segment)
        end
    end

    # Then print report text
    open(reportpath, "w") do reportfile
        for state in states
            println(reportfile, report_text(state), '\n')
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 4
        error("Usage: julia report.jl conspaths matpaths accpath reportpath")
    end
    conspaths = split(ARGS[1], ',')
    matpaths = split(ARGS[2], ',')
    length(conspaths) == length(matpaths) || error("must have same number of cons and matpaths")
    main(conspaths, matpaths, ARGS[3], ARGS[4])
end
