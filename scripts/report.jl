# TODO: Check for extra nt at ends, by aligning to a selection of conserved ends

module t

using FASTX
using BioSequences
using BioAlignments
using CodecZlib
using ErrorTypes
using Printf
using Serialization

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

const ALN_MODEL = AffineGapScoreModel(EDNAFULL, gap_open=-25, gap_extend=-2)

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

#=
1. Overlap alignment of segment to ORF
2. Extract 
=#

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
                # Headers look like "HEADER_HA"
                p = findlast(isequal('_'), line)
                p === nothing && error("Found header \"$(line)\", expected pattern \"HEADER_SEGMENT\"")
                segment = Segment(line[p+1:end])
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

struct ORF
    var::ProteinVariant
    seqs::Union{Tuple{LongDNASeq}, NTuple{2, LongDNASeq}}
end

# TODO: Too inefficient to load this again and again and again...right?
# maybe instead do ALL the reports in one go, instead of individual reports
# then gathering them.
function load_orfs(refdir::AbstractString, segment::Segment,
accession::AbstractString)::Union{Tuple{ORF}, NTuple{2, ORF}}
    record = FASTA.Record()
    fastapath = joinpath(refdir, string(segment) * ".fna")
    seq = open(FASTA.Reader, fastapath) do reader
        while !eof(reader)
            read!(reader, record)
            if FASTA.identifier(record) == accession
                return FASTA.sequence(LongDNASeq, record)
            end
        end
        error("Accession $accession not found in $fastapath")
    end
    orfpath = joinpath(refdir, string(segment) * ".orfs.jls")
    open(orfpath) do io
        for (acc, v) in deserialize(io)
            if acc == accession
                result = ORF[]
                for protein in v
                    protseqs = Tuple([seq[orf] for orf in protein[2]])
                    push!(result, ORF(ProteinVariant(protein[1]), protseqs))
                end
                return Tuple(result)
            end
        end
        error("Accession $accession not found in $orfpath")
    end
end

#function foo(orf::ORF, segment::LongDNASeq)
#    aln = pairalign(OverlapAlignment(), segment, orf.seq, ALN_MODEL).aln



struct Assembly
    segment::Segment
    insignificant::BitVector
    seq::LongDNASeq
    orfs::Union{Tuple{ORF}, NTuple{2, ORF}}
end

# We only have this to reconstruct the record - including the lowercase letters.
function FASTA.Record(x::Assembly, basename::AbstractString)
    data = UInt8['>'; codeunits(basename); '_'; codeunits(string(x.segment)); '\n']
    seqdata = take!((io = IOBuffer(); print(io, x.seq); io))
    @assert length(seqdata) == length(x.insignificant)
    for (i, is_insignificant) in enumerate(x.insignificant)
        is_insignificant && (seqdata[i] += 0x20)
    end
    FASTA.Record(append!(data, seqdata))
end

function load_assembly(assemblypath::AbstractString, orfdir::AbstractString)::Option{Vector{Assembly}}
    records = open(FASTA.Reader, assemblypath) do reader
        map(reader) do record
            v = rsplit(FASTA.identifier(record), '_', limit=2)
            length(v) == 2 || error("Found header \"$(id)\", expected pattern \"HEADER_SEGMENT\"")
            accession = first(v)
            segment = Segment(last(v))
            seq = FASTA.sequence(LongDNASeq, record)
            @assert length(record.sequence) == length(seq)
            insignificant = BitVector([in(i, UInt8('a'):UInt8('z')) for i in @view record.data[record.sequence]])
            orfs = load_orfs(orfdir, segment, accession)
            Assembly(segment, insignificant, seq, orfs)
        end
    end
    # Note: MUST be sorted!
    isempty(records) ? none : Thing(sort!(records, by=x -> x.segment))
end

####################################
function report_text(assemblies::Vector{Assembly}, matpath::AbstractString)::String
    # Get the stuff from the input files
    names_depths = expect(get_depths(matpath), "Empty matrix file: $matpath")
    if map(first, names_depths) != [a.segment for a in assemblies]
        error("Differing segments in matrix path compared to assembly: $matpath")
    end

    lines = String[]
    for segment in instances(Segment)
        i = findfirst(x -> first(x) == segment, names_depths)
        data::Option{Tuple{Vector{UInt32}, Assembly}} = if i === nothing
            none
        else
            Thing((names_depths[i][2], assemblies[i]))
        end
        append!(lines, report_text(segment, data))
    end
    join(lines, '\n')
end

function report_text(segment::Segment, data::Option{Tuple{Vector{UInt32}, Assembly}})::Vector{String}
    # Check segment length and return early if way too short
    is_none(data) && return ["$segment:\n\tERROR: Missing segment"]
    depths, assembly = unwrap(data)

    minlen = min(length(assembly.seq), length(depths))
    if minlen < 2 * TERMINAL + 1
        return ["$segment:\n\tERROR: Sequence or mat.gz length: $minlen"]
    end

    # Header: HA: depth 4.32e+03 coverage 1.000 
    mean_depth = sum(UInt, depths) / length(depths)
    coverage = count(!iszero, depths) / length(depths)
    result = ["$(rpad(string(segment) * ':', 4)) depth $(@sprintf "%.2e" mean_depth) coverage $(@sprintf "%.3f" coverage)"]

    # insignificant bases
    n_insignificant = count(assembly.insignificant)
    if !iszero(n_insignificant)
        push!(result, "\t       $(n_insignificant) bases insignificantly basecalled")
    end

    # Check for low depths
    lowdepth = count(x -> x < 25, @view depths[TERMINAL + 1: end - TERMINAL])
    if !iszero(lowdepth)
        push!(result, "\t       $lowdepth central bases with depth < 25")
    end

    # Check for ambiguous bases
    amb = count(isambiguous, assembly.seq)
    if !iszero(amb)
        push!(result, "\t       $amb ambiguous bases in consensus sequence")
    end
    
    result
end

function main(outdir::AbstractString, basename::AbstractString,
    assemblypath::AbstractString, matpath::AbstractString, orfdir::AbstractString
)
    assemblies = expect(load_assembly(assemblypath, orfdir), "Empty kma assembly file: $assemblypath")

    # Make report
    open(joinpath(outdir, "report.txt"), "w") do reportfile
        println(reportfile, report_text(assemblies, matpath))
    end
    
    # Create consensus files
    open(FASTA.Writer, joinpath(outdir, "consensus.fna")) do writer
        for assembly in assemblies
            write(writer, FASTA.Record(assembly, basename))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 4
        error("Usage: julia report.jl outdir basename assemblypath matrixpath")
    end
    main(ARGS...)
end

end # module