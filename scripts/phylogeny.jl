# Directory structure of ref/trees is:

# . refset [e.g. human / swine / avian]
# ├── segment1 [segment name: e.g. HA]
# │   ├── subtype1.fna [arbitrary: e.g. H1sw.fna]
# │   └── subtype2.fna [arbitrary: e.g. H2.fna] 
# └── segment2 [e.g. NA]
#     ├── subtype1.fna [arbitrary: e.g. N1.fna]
#     └── subtype2.fna [arbitrary: e.g. N2sw.fna]

# Each segment is processed independently.
# * First, the closest subtype is found by roughly comparing against the sequences
# in each subtype.
# * Then, the segment is added to the subtype fasta to obtain the "cat"
# * Then, the cat is aligned and trimmed
# * Then, the trimmed is run through phylogeny

using FASTX
using InfluenzaCore: Segment
using BioSequences: LongDNASeq, DNAMer, each
using ErrorTypes
using MinHash

function main(
    outdir::String,
    refoutdir::String,
    refdir::String,
    consensusdir::String
)
    # Get segments and possible subtypes
    refs = load_refs(refdir)

    # Figure out which sequences have passed quality control
    segments = Tuple(keys(refs))
    consensus = load_consensus(segments, consensusdir)

    # Get best subtype for each
    best_subtypes = get_subtype(segments, consensus, refs)

    # Make sure guide trees for each subtype has been calculated
    subtypes = Dict(seg => collect(keys(subt)) for (seg, subt) in refs)
    create_guide_trees(subtypes, refoutdir)


    # Align consensus sequences to reference alignment
    # Run IQ-TREE
    nothing
end

struct Basename
    name::String
end

struct FluType
    name::String
end

for T in (Basename, FluType)
    @eval begin
        name(x::$T) = x.name
        Base.hash(x::$T, h::UInt) = hash(name(x), h)
        Base.:(==)(x::$T, y::$T) = name(x) == name(y)
    end
end

const N_SEGMENTS = length(instances(Segment))
const SegmentTuple{T} = NTuple{N_SEGMENTS, T}
const MINHASHER = MinHasher(250)
const K = 7

"Open each ref, and get the flu subtypes, and min hash sketches for each of these"
function load_refs(refseqdir::String
)::Dict{Segment, Dict{FluType, Vector{Tuple{LongDNASeq, MinHashSketch}}}}
    record = FASTA.Record()
    segments = map(readdir(refseqdir)) do subdir
        @unwrap_or parse(Segment, subdir) error("Could not parse directory as segment: $subdir")
    end
    result = Dict{Segment, Dict{FluType, Vector{Tuple{LongDNASeq, MinHashSketch}}}}()
    for segment in segments
        result[segment] = valtype(result)()
        subtypes = map(readdir(joinpath(refseqdir, string(segment)))) do filename
            @assert endswith(filename, ".fna")
            FluType(filename[1:end-4])
        end
        for subtype in subtypes
            result[segment][subtype] = valtype(valtype(result))()
            open(FASTA.Reader, joinpath(refseqdir, string(segment), name(subtype) * ".fna")) do reader
                while !eof(reader)
                    read!(reader, record)
                    empty!(MINHASHER)
                    seq = FASTA.sequence(LongDNASeq, record)
                    update!(MINHASHER, each(DNAMer{K}, seq))
                    push!(result[segment][subtype], (seq, MinHashSketch(MINHASHER)))
                end
            end
        end
    end
    return result
end

"Load all sequences. The vector in the output values are in the same order as
the input segments vector"
function load_consensus(segments::NTuple{N, Segment}, consensusdir::String
)::Dict{Basename, NTuple{N, Option{Tuple{LongDNASeq, MinHashSketch}}}} where N
    isdir(consensusdir) || error("No such directory: $consensusdir")
    result = Dict{Basename, NTuple{N, Option{Tuple{LongDNASeq, MinHashSketch}}}}()
    record = FASTA.Record()
    buffer = fill(none(Tuple{LongDNASeq, MinHashSketch}), N)
    for dir in readdir(consensusdir)
        fill!(buffer, none(Tuple{LongDNASeq, MinHashSketch}))
        basename = Basename(dir)
        filename = joinpath(consensusdir, dir, "curated.fna")
        if isfile(filename)
            open(FASTA.Reader, filename) do reader
                while !eof(reader)
                    read!(reader, record)
                    segment = unwrap(parse(Segment, last(rsplit(FASTA.header(record), '_', limit=2))))
                    index = findfirst(isequal(segment), segments)
                    index === nothing && continue
                    seq = FASTA.sequence(LongDNASeq, record)
                    empty!(MINHASHER)
                    update!(MINHASHER, each(DNAMer{K}, seq))
                    buffer[index] = some((seq, MinHashSketch(MINHASHER)))
                end
            end
        end
        result[basename] = NTuple{N}(buffer)
    end
    return result
end

"Get best subtype for each"

function get_subtype(
    segments::NTuple{N, Segment},
    consensus::Dict{Basename, NTuple{N, Option{Tuple{LongDNASeq, MinHashSketch}}}},
    refs::Dict{Segment, Dict{FluType, Vector{Tuple{LongDNASeq, MinHashSketch}}}},
)::Dict{Basename, NTuple{N, Option{FluType}}} where N
    result = Dict{Basename, NTuple{N, Option{FluType}}}()
    buffer = fill(none(FluType), N)
    for (basename, maybes) in consensus
        fill!(buffer, none(FluType))
        for (i, segment, maybe_pair) in zip(1:N, segments, maybes)
            (seq, seq_minhash) = @unwrap_or maybe_pair continue
            best_score = 0
            best_subtype = nothing
            for (subtype, seqpairs) in refs[segment], (seq, refseq_minhash) in seqpairs
                score = intersectionlength(seq_minhash, refseq_minhash)
                if score > best_score
                    best_score = score
                    best_subtype = subtype
                end
            end
            if best_subtype !== nothing
                buffer[i] = some(best_subtype)
            end
        end
        result[basename] = NTuple{N}(buffer)
    end
    return result
end

function create_guide_trees(
    subtypes::Dict{Segment, Vector{FluType}},
    refs::Dict{Segment, Dict{FluType, Vector{Tuple{LongDNASeq, MinHashSketch}}}},
    refoutdir::String
)
    for (segment, types) in subtypes, type in types
        seqs = map(first, refs[segment][flutype])
        outdir = joinpath(refoutdir, segment)
        create_guide_tree(segment, type, seqs, outdir)
    end
end

function create_guide_tree(
    segment::Segment, flutype::FluType,
    seqs::Vector{LongDNASeq}, outdir::String
)
    trimpath = joinpath("outdir", string(segment) * ".aln.trim.fna")
    treepath = joinpath("outdir", string(segment) * ".treefile")
    ispath(treepath) && return nothing
    
    # Align if alignment doesn't exist



    # Trim if trim doesn't exist
    # Create guidetree if it doesn't exist
end