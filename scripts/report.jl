#=
What it does
* Read in references from aln/cat.fna
* Get the ORFs from ref/XXXXXX


=#

using FASTX
using BioSequences
using BioAlignments
using CodecZlib
using ErrorTypes
using Printf

# Number of nucleotides from each end to consider the "ends" of the sequence
# in which we have more lax requirements for depth
const TERMINAL = 25

"""
    Segment(::String) -> Segment

Create a representaiton of one of the eight Influenza A genome segments.
"""
@enum Segment::UInt8 begin
    PB1
    PB2
    PA
    HA
    NP
    NA
    MP
    NS
end

const _STR_SEGMENT = Dict(map(i -> string(i)=>i, instances(Segment)))
function Segment(s::AbstractString)
    smt = get(_STR_SEGMENT, s, nothing)
    smt === nothing && error("Unknown segment: \"$s\"")
    smt
end

@enum ProteinVariant::UInt8 begin
    pb2
    pb1
    pb1fa
    pa
    pax
    ha
    np
    na
    m1
    m2
    ns1
    nep
end

const _CRITICAL = [true, true, false, true, false, 
    true, true, true, true, true, true, true]
is_critical(x::ProteinVariant) = @inbounds _CRITICAL[reinterpret(UInt8, x) + 0x01]

struct ORF
    variant::ProteinVariant
    segment::LongDNASeq
    exons::Vector{UnitRange{UInt16}}
end

const ALN_MODEL = AffineGapScoreModel(EDNAFULL, gap_open=-13, gap_extend=-2)

function Base.join(v::AbstractVector{T}) where T <: LongSequence
    ln = sum(length, v)
    seq = T(ln)
    offset = 0
    @inbounds for s in v
        seq[offset+1:offset+length(s)] = s
        offset += length(s)
    end
    seq
end

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

"Given a path to a .mat.gz file, return a vector of (segment, depths)"
function get_depths(matpath::AbstractString)::Option{Vector{Tuple{Segment, Vector{UInt32}}}}
    open(matpath) do io
        result = Vector{Tuple{Segment, Vector{UInt32}}}()
        segment = nothing
        depths = UInt32[]
        fields = Vector{SubString{String}}(undef, 7)
        linedepths = Vector{UInt32}(undef, 6)
        for line in Iterators.map(strip, eachline(GzipDecompressorStream(io)))
            if isempty(line)
                if !isempty(depths)
                    push!(result, (segment, depths))
                    depths = UInt32[]
                end
                continue
            end
            if startswith(line, '#')
                segment = Segment(line[2:end])
                continue
            end
            split!(fields, line, UInt8('\t'))
            @inbounds for i in 1:6
                linedepths[i] = parse(UInt32, fields[i+1])
            end
            # Skip lines with majority as deletion
            argmax(linedepths) == 6 && continue
            total_depth = sum(linedepths)
            # We don't include N in depths!
            depth = total_depth - linedepths[5]
            push!(depths, depth)
        end
        isempty(depths) || push!(result, (segment, depths))
        isempty(result) && return none
        return Thing(sort!(result))
    end
end

"Loads sequences, returns a vector like: [(HA, TAGTAGTCGAT, Thing(2)), ... ], last is n_insignificant
inside the terminals"
function load_consensus(consensuspath::AbstractString)::Option{Vector{Tuple{Segment, LongDNASeq, Option{UInt32}}}}
    records = open(FASTA.Reader, consensuspath) do reader
        map(reader) do record
            n_insignificant = if length(record.sequence) < (2 * TERMINAL + 1)
                none
            else
                n = count(@view record.data[record.sequence][TERMINAL + 1 : end-TERMINAL]) do byte
                    in(byte, UInt8('a'):UInt('z'))
                end
                Thing(UInt32(n))
            end
            (Segment(FASTA.identifier(record)),  FASTA.sequence(LongDNASeq, record), n_insignificant)
        end
    end
    isempty(records) ? none : Thing(sort!(records))
end

####################################
function report_text(names_seqs::Vector{Tuple{Segment, LongDNASeq, Option{UInt32}}}, matpath::AbstractString)::String
    # Get the stuff from the input files
    names_depths = expect(get_depths(matpath), "Empty matrix file: $matpath")
    if map(first, names_depths) != map(first, names_seqs)
        error("Differing segments in files $conspath, $matpath")
    end

    lines = String[]
    for ((segment, depths), (_, seq, n_insignificant)) in zip(names_depths, names_seqs)
        append!(lines, report_text(segment, depths, seq, n_insignificant))
    end
    join(lines, '\n')
end

function report_text(segment::Segment, depths::Vector{<:Unsigned}, seq::LongDNASeq,
    n_insignificant::Option{<:Unsigned}
)::Vector{String}
    #length(seq) == length(depths) || error("Unequal length for seq and depths for $segment")
    
    # Check segment length and return early if way too short
    minlen = min(length(seq), length(depths))
    if minlen < 2 * TERMINAL + 1
        return ["$segment:\n\tERROR: Sequence or mat.gz length: $minlen"]
    end

    # Header: HA: depth 4.32e+03 coverage 1.000 
    mean_depth = sum(UInt, depths) / length(depths)
    coverage = count(!iszero, depths) / length(depths)
    result = ["$(rpad(string(segment) * ':', 4)) depth $(@sprintf "%.2e" mean_depth) coverage $(@sprintf "%.3f" coverage)"]

    # insignificant bases
    if !iszero(unwrap(n_insignificant))
        push!(result, "\t       $(unwrap(n_insignificant)) bases insignificantly basecalled")
    end

    # Check for low depths
    lowdepth = count(x -> x < 25, @view depths[TERMINAL + 1: end - TERMINAL])
    if !iszero(lowdepth)
        push!(result, "\t       $lowdepth central bases with depth < 25")
    end

    # Check for ambiguous bases
    amb = count(isambiguous, seq[TERMINAL + 1: end - TERMINAL])
    if !iszero(amb)
        push!(result, "\t       $amb central ambiguous bases in consensus sequence")
    end
    
    result
end

function main(outdir::AbstractString, basename::AbstractString,
    assemblypath::AbstractString, matpath::AbstractString
)
    names_seqs = expect(load_consensus(assemblypath), "Empty kma assembly file: $assemblypath")

    # Make report
    open(joinpath(outdir, "report.txt"), "w") do reportfile
        println(reportfile, report_text(names_seqs, matpath))
    end
    
    # Create consensus files
    open(FASTA.Writer, joinpath(outdir, "consensus.fna")) do writer
        for (segment, seq) in names_seqs
            header = basename * '_' * string(segment)
            write(writer, FASTA.Record(header, seq))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 4
        error("Usage: julia report.jl outdir basename assemblypath matrixpath")
    end
    main(ARGS...)
end
