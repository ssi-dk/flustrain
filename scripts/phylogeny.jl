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
using BioSequences: LongDNASeq, DNAMer, each, NucleicAcidAlphabet, BioSequence, fwmer
using ErrorTypes
using MinHash

function main(
    outdir::String,
    refoutdir::String, # e.g. refout/human
    refdir::String, # e.g. ref/trees/human
    consensusdir::String
)
    # Get segments and possible subtypes
    refpaths = get_refpaths(refdir)
    ref_sketches = load_ref_sketches(refpaths)

    # Figure out which sequences have passed quality control
    segments = Tuple(keys(ref_sketches))
    cons_seqs, cons_sketches = load_consensus(segments, consensusdir)

    # Make paths
    foreach(bn -> mkpath(joinpath(outdir, name(bn))), keys(cons_seqs))

    # Get best subtype for each
    best_subtypes = get_subtype(segments, cons_sketches, ref_sketches)

    # Write out subtypes
    write_subtypes(segments, best_subtypes, outdir)
    
    # Make sure guide trees for each subtype has been calculated
    subtypes = Dict(seg => collect(keys(subt)) for (seg, subt) in ref_sketches)
    create_guide_trees(subtypes, refdir, refoutdir)

    # Align consensus sequences to reference alignment
    align_consensus(segments, cons_seqs, best_subtypes, refoutdir, outdir)

    # Run IQ-TREE
    iq_tree(segments, best_subtypes, refoutdir, outdir)
    return nothing
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

function MinHash.update!(mh::MinHasher, seq::BioSequence{<:NucleicAcidAlphabet})
    update!(mh, (fwmer(mer) for mer in each(DNAMer{K}, seq)))
end
MinHash.update!(mh::MinHasher, rec::FASTA.Record) = update!(mh, FASTA.sequence(LongDNASeq, rec))

function get_refpaths(refseqdir::String)::Dict{Segment, Dict{FluType, String}}
    result = Dict{Segment, Dict{FluType, String}}()
    segments = map(readdir(refseqdir)) do subdir
        @unwrap_or parse(Segment, subdir) error("Could not parse directory as segment: $subdir")
    end
    for segment in segments
        result[segment] = valtype(result)()
        dirpath = joinpath(refseqdir, string(segment))
        for filename in readdir(dirpath)
            endswith(filename, ".fna") || error("Filename should end with .fna: $filename")
            subtype = FluType(filename[1:end-4])
            path = joinpath(dirpath, filename)
            result[segment][subtype] = path
        end
    end
    return result
end

"Open each ref, and get the flu subtypes, and min hash sketches for each of these"
function load_ref_sketches(refpaths::Dict{Segment, Dict{FluType, String}}
)::Dict{Segment, Dict{FluType, Vector{MinHashSketch}}}
    record = FASTA.Record()
    result = Dict(s => Dict(f => MinHashSketch[] for f in keys(fs)) for (s, fs) in refpaths)
    for (segment, subtypedir) in refpaths, (subtype, path) in subtypedir
        open(FASTA.Reader, path) do reader
            while !eof(reader)
                read!(reader, record)
                empty!(MINHASHER)
                update!(MINHASHER, record)
                push!(result[segment][subtype], MinHashSketch(MINHASHER))
            end
        end
    end
    return result
end

"Load all sequences. The vector in the output values are in the same order as
the input segments vector"
function load_consensus(segments::NTuple{N, Segment}, consensusdir::String
)::Tuple{
    Dict{Basename, NTuple{N, Option{LongDNASeq}}} where N,
    Dict{Basename, NTuple{N, Option{MinHashSketch}}} where N,
} where N
    isdir(consensusdir) || error("No such directory: $consensusdir")
    seqs = Dict{Basename, NTuple{N, Option{LongDNASeq}}}()
    minhashes = Dict{Basename, NTuple{N, Option{MinHashSketch}}}()
    record = FASTA.Record()
    seqbuffer = fill(none(LongDNASeq), N)
    sketchbuffer = fill(none(MinHashSketch), N)
    for dir in readdir(consensusdir)
        fill!(seqbuffer, none(LongDNASeq))
        fill!(sketchbuffer, none(MinHashSketch))
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
                    update!(MINHASHER, seq)
                    seqbuffer[index] = some(seq)
                    sketchbuffer[index] = some(MinHashSketch(MINHASHER))
                end
            end
        end
        seqs[basename] = NTuple{N}(seqbuffer)
        minhashes[basename] = NTuple{N}(sketchbuffer)
    end
    return (seqs, minhashes)
end

"Get best subtype for each"
function get_subtype(
    segments::NTuple{N, Segment},
    cons_sketches::Dict{Basename, NTuple{N, Option{MinHashSketch}}},
    refsketches::Dict{Segment, Dict{FluType, Vector{MinHashSketch}}},
)::Dict{Basename, NTuple{N, Option{FluType}}} where N
    result = Dict{Basename, NTuple{N, Option{FluType}}}()
    buffer = fill(none(FluType), N)
    for (basename, maybe_sketches) in cons_sketches
        fill!(buffer, none(FluType))
        for (i, segment, maybe_sketch) in zip(1:N, segments, maybe_sketches)
            seq_minhash = @unwrap_or maybe_sketch continue
            best_score = 0
            best_subtype = nothing
            for (subtype, minhashes) in refsketches[segment], ref_minhash in minhashes
                score = intersectionlength(seq_minhash, ref_minhash)
                if score > best_score
                    best_score = score
                    best_subtype = subtype
                end
            end
            if best_subtype !== nothing
                buffer[i] = some(best_subtype)
            else
                @warn "Basename $(name(basename)) segment $(string(segment)) does not look like any subtype"
            end
        end
        result[basename] = NTuple{N}(buffer)
    end
    return result
end

function write_subtypes(
    segments::NTuple{N, Segment},
    best_subtypes::Dict{Basename, NTuple{N, Option{FluType}}},
    outdir::String
) where N
    for (basename, maybe_types) in best_subtypes
        path = joinpath(outdir, name(basename), "subtypes.txt")
        open(path, "w") do file
            for (segment, maybe_type) in zip(segments, maybe_types)
                subtype = @unwrap_or maybe_type continue
                println(file, string(segment), '\t', name(subtype))
            end
        end
    end
end

function create_guide_trees(
    subtypes::Dict{Segment, Vector{FluType}},
    refseqdir::String,
    refoutdir::String
)
    for (segment, types) in subtypes, type in types
        refseqpath = joinpath(refseqdir, string(segment), name(type) * ".fna")
        outdir = joinpath(refoutdir, string(segment))
        create_guide_tree(type, refseqpath, outdir)
    end
end

function create_guide_tree(flutype::FluType, refseqpath::String, outdir::String)
    trimpath = joinpath(outdir, name(flutype) * ".aln.trim.fna")
    if !ispath(trimpath)
        mkpath(outdir)
        alnpath = joinpath(outdir, name(flutype) * ".aln.fna")
        run(pipeline(`mafft --thread 2 $refseqpath`, stdout=alnpath, stderr=devnull))
        cmd = `trimal -gt 0.9 -cons 60 -keepheader -keepseqs -in $alnpath -fasta -out $trimpath`
        run(cmd)
    end

    # IQTREE if the guide tree doesn't exist
    treefilepath = joinpath(outdir, name(flutype) * ".treefile")
    if !ispath(treefilepath)
        prepath = joinpath(outdir, name(flutype))
        cmd = `iqtree -T 4 --quiet -s $trimpath -pre $prepath -m HKY+G2`
        run(cmd)
    end
end

function align_consensus(
    segments::NTuple{N, Segment},
    consensus::Dict{Basename, NTuple{N, Option{LongDNASeq}}},
    best_subtypes::Dict{Basename, NTuple{N, Option{FluType}}},
    refoutdir::String,
    outdir::String
) where N
    tasks = Task[]
    for (basename, seqtuple) in consensus
        for (segment, maybe_type, maybe_seq) in zip(segments, best_subtypes[basename], consensus[basename])
            seq = @unwrap_or maybe_seq continue
            subtype = @unwrap_or maybe_type continue
            basedir = joinpath(outdir, name(basename))
            treealnpath = joinpath(refoutdir, string(segment), name(subtype) * ".aln.trim.fna")
            @assert ispath(treealnpath)
            push!(tasks, Threads.@spawn align_consensus(segment, seq, basedir, treealnpath))
        end
    end
    foreach(wait, tasks)
    return nothing
end

function align_consensus(
    segment::Segment, seq::LongDNASeq, basedir::String, treealnpath::String
)
    # Write out each present consensus sequence.
    seqpath = joinpath(basedir, string(segment) * ".fna")
    open(seqpath, "w") do file
        println(file, ">$(string(segment))\n", seq)
    end

    # Align
    alnpath = joinpath(basedir, string(segment) * ".aln.fna")
    ispath(alnpath) || run(pipeline(`mafft --add $seqpath --keeplength $treealnpath`, stdout=alnpath))
    return nothing
end

function iq_tree(
    segments::NTuple{N, Segment},
    best_subtypes::Dict{Basename, NTuple{N, Option{FluType}}}, # use this to check which aln are present
    refoutdir::String,
    outdir::String
) where N
    tasks = Task[]
    for (basename, maybe_types) in best_subtypes, (segment, maybe_type) in zip(segments, maybe_types)
        subtype = @unwrap_or maybe_type continue
        alnpath = joinpath(outdir, name(basename), string(segment) * ".aln.fna")
        prepath = joinpath(outdir, name(basename), string(segment))
        target_path = prepath * ".iqtree"
        guidepath = joinpath(refoutdir, string(segment), name(subtype) * ".treefile")
        @assert isfile(guidepath)
        @assert isfile(alnpath)
        cmd = `iqtree -T 1 -g $guidepath --quiet -s $alnpath -pre $prepath -m HKY+G2`
        isfile(target_path) || push!(tasks, Threads.@spawn run(cmd))
    end
    foreach(wait, tasks)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 4
        println("Usage: julia phylogeny.jl outdir refoutdir refdir consensusdir")
        exit(1)
    else
        main(ARGS...)
    end
end