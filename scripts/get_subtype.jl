using FASTX
using InfluenzaCore: Segment
using BioSequences: LongDNASeq, DNAMer, each
using ErrorTypes
using MinHash

struct SubType
    name::String
end
name(s::SubType) = s.name

struct Basename
    name::String
end
name(s::Basename) = s.name

const MINHASHER = MinHasher(250)
const RECORD = FASTA.Record()
const K = 7

function main(refseqdir::String, consdir::String, outdir::String)
    # Load the tree reference sequences
    refs = load_refs(refseqdir)

    # Load all curated sequences
    consensus = load_consensus(consdir, Set(keys(refs)))

    # Do some rough measure to establish the best subtype for each consensus seq
    best = get_subtype(consensus, refs)

    # Write out sequences and best
    write_results(consensus, best, outdir)
    return nothing
end

function load_refs(refseqdir::String)::Dict{Segment, Dict{SubType, Vector{MinHashSketch}}}
    segments = map(readdir(refseqdir)) do subdir
        @unwrap_or parse(Segment, subdir) error("Could not parse directory as segment: $subdir")
    end
    result = Dict{Segment, Dict{SubType, Vector{MinHashSketch}}}()
    for segment in segments
        result[segment] = Dict{SubType, Vector{MinHashSketch}}()
        subtypes = map(readdir(joinpath(refseqdir, string(segment)))) do filename
            @assert endswith(filename, ".fna")
            SubType(filename[1:end-4])
        end
        for subtype in subtypes
            result[segment][subtype] = MinHashSketch[]
            open(FASTA.Reader, joinpath(refseqdir, string(segment), name(subtype) * ".fna")) do reader
                while !eof(reader)
                    read!(reader, RECORD)
                    empty!(MINHASHER)
                    update!(MINHASHER, each(DNAMer{K}, FASTA.sequence(LongDNASeq, RECORD)))
                    push!(result[segment][subtype], MinHashSketch(MINHASHER))
                end
            end
        end
    end
    return result
end

function load_consensus(consensusdir::String, segments::Set{Segment}
)::Dict{Basename, Dict{Segment, Option{Tuple{LongDNASeq, MinHashSketch}}}}
    result = Dict{Basename, Dict{Segment, Option{Tuple{LongDNASeq, MinHashSketch}}}}()
    for basename in map(Basename, readdir(consensusdir))
        result[basename] = Dict(segment => none(Tuple{LongDNASeq, MinHashSketch}) for segment in segments)
        subdir = joinpath(consensusdir, name(basename))
        filename = joinpath(subdir, "curated.fna")
        isfile(filename) || continue
        open(FASTA.Reader, filename) do reader
            while !eof(reader)
                read!(reader, RECORD)
                segment = unwrap(parse(Segment, rsplit(FASTA.header(RECORD), '_', limit=2)[2]))
                in(segment, segments) || continue
                empty!(MINHASHER)
                seq = FASTA.sequence(LongDNASeq, RECORD)
                update!(MINHASHER, each(DNAMer{K}, seq))
                result[basename][segment] = some((seq, MinHashSketch(MINHASHER)))
            end
        end
    end
    return result
end

function get_subtype(
    consensus::Dict{Basename, Dict{Segment, Option{Tuple{LongDNASeq, MinHashSketch}}}},
    refs::Dict{Segment, Dict{SubType, Vector{MinHashSketch}}}
)::Dict{Basename, Dict{Segment, Option{SubType}}}
    result = Dict{Basename, Dict{Segment, Option{SubType}}}()
    for (basename, seqdict) in consensus
        result[basename] = Dict(s => none(SubType) for s in keys(seqdict))
        for (segment, maybe_seq) in seqdict
            (seq, seq_minhash) = @unwrap_or maybe_seq continue
            best_score = 0
            best_subtype = nothing
            for (subtype, refseq_minhashes) in refs[segment]
                for refseq_minhash in refseq_minhashes
                    score = intersectionlength(seq_minhash, refseq_minhash)
                    if score > best_score
                        best_score = score
                        best_subtype = subtype
                    end
                end
            end
            best_subtype === nothing && error("$(string(segment)) of $basename doesn't look like any subtype")
            result[basename][segment] = some(best_subtype)
        end
    end
    return result
end

function write_results(
    consensus::Dict{Basename, Dict{Segment, Option{Tuple{LongDNASeq, MinHashSketch}}}},
    best::Dict{Basename, Dict{Segment, Option{SubType}}},
    outdir::String
)::Nothing
    for (basename, seqdict) in consensus
        for (segment, maybe_seq) in seqdict
            (seq, seq_minhash) = @unwrap_or maybe_seq continue
            open(joinpath(outdir, name(basename), string(segment) * ".fna"), "w") do file
                println(file, '>', name(basename), '_', string(segment))
                println(file, seq)
            end
            open(joinpath(outdir, name(basename), string(segment) * "_subtype.txt"), "w") do file
                println(file, name(unwrap(best[basename][segment])))
            end
        end
    end
    return nothing
end

if length(ARGS) != 3
    error("Usage: julia get_subtype.jl refseqdir consdir outdir")
end
main(ARGS...)
s