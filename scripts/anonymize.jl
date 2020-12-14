# This script removes all reads that do not map to a given reference

using CodecZlib
using FASTX
using SHA
using Comonicon

const ENDINGS = ("fastq.gz", ".fq.gz")
const PATTERN = r"^([^_]+_S\d+)_L001_R(1|2)_001.fastq.gz$"

struct Hash # 128 bits is enough
    a::UInt64
    b::UInt64
end

function Hash(v::Union{AbstractString, AbstractArray{UInt8}})
    sh = sha256(v)
    GC.@preserve sh begin
        p = Ptr{UInt64}(pointer(sh))    
        h = Hash(unsafe_load(p, 1), unsafe_load(p, 2))
    end
    return h
end

function Hash(rec::FASTQ.Record)
    isempty(rec.identifier) && throw(ArgumentError("Read empty"))
    Hash(view(rec.data, rec.identifier))
end

##
function parse_filename(s::AbstractString)
    m = match(PATTERN, s)
    m === nothing && throw(ArgumentError("Invalid file name: $s"))
    head = m[1]
    fw = m[2] == "1"
    return head, fw
end

function get_fastqs(dir::AbstractString)::Vector{Tuple{String, String}}
    files = readdir(dir)
    
    # Remove all that are not fastq files
    sort!(filter!(x -> any(endswith(x, i) for i in ENDINGS), files))
    d = Dict()

    for f in files
        head, isfw = parse_filename(f)
        v = get(d, head, nothing)
        if v === nothing
            d[head] = [f]
        elseif length(v) > 1
            error("Sample $head present in >2 copies")
        else
            push!(v, f)
        end
    end

    single = findfirst(x -> length(x) == 1, collect(values(d)))
    single === nothing || error("Sole sample: $single")
    return [Tuple(sort(v)) for v in values(d)]
end

"Parse the out.frag.gz file and extract a set of mapping readnames"
function build_present_set(path::AbstractString, minscore::Int)::Set{Hash}
    open(path) do file
        io = GzipDecompressorStream(file)
        s = sizehint!(Set{Hash}(), 100000)
        for (lineno, line) in enumerate(eachline(io))
            stripped = strip(line)
            isempty(stripped) && continue
            fields = split(stripped, '\t')
            if (length(fields) != 7) & (length(fields) != 9)
                throw(ArgumentError("Line $lineno of file $path does not contain 7 or 9 fields"))
            end
            score = parse(Int, fields[3])
            identifier = first(split(last(fields)))
            score > minscore && push!(s, Hash(identifier))
        end
        return s
    end
end

function filter_fastq(s::Set{Hash}, inp::AbstractString, out::AbstractString)
    reader = FASTQ.Reader(GzipDecompressorStream(open(inp)))
    writer = FASTQ.Writer(GzipCompressorStream(open(out, "w")))
    rec = FASTQ.Record()
    while !eof(reader)
        read!(reader, rec)
        if in(Hash(rec), s)
            write(writer, rec)
        end
    end
    close(reader)
    close(writer)
end

##
function kma_map(dir::AbstractString, fw::AbstractString, rv::AbstractString, db::AbstractString)
    # Memory map input (we assume small files), and do not create consensus seqs
    run(`kma -ipe $fw $rv -o $dir/kmaout -mmap -nc -t_db $db`)
end

function filter_fastq(dstdir, srcdir, fw, rv, db, minscore)
    fwin = joinpath(srcdir, fw)
    fwout = joinpath(dstdir, fw)
    rvin = joinpath(srcdir, rv)
    rvout = joinpath(dstdir, rv)

    hashset = mktempdir(".") do dir
        kma_map(dir, fwin, rvin, db)
        hashset = build_present_set(joinpath(dir, "kmaout.frag.gz"), minscore)
    end

    filter_fastq(hashset, fwin, fwout)
    filter_fastq(hashset, rvin, rvout)
end

function _anonymize(dstdir, srcdir, kmadb, minscore)
    mkdir(dstdir)
    pairs = get_fastqs(srcdir)
    Threads.@threads for (fw, rv) in pairs
        filter_fastq(dstdir, srcdir, fw, rv, kmadb, minscore)
    end
end

"""
Anonymize reads not mapping to KMA index.

This script assumes
* Short, high-quality, paired-end Illumina-like reads.
* Reads and index small enough to align in memory
* gzipped FASTQ filenames match regex `$PATTERN`

# Arguments

- `dstdir`: Destination directory to create
- `srcdir`: Where to find the gzipped FASTQ files
- `kmadb`: Name of KMA index, without file extention

# Options

- `--minscore <arg>`: Minimum alignment score to keep [100]
"""
if abspath(PROGRAM_FILE) == @__FILE__
    @main function anonymize(dstdir, srcdir, kmadb, minscore::Int=100)
        _anonymize(dstdir, srcdir, kmadb, minscore)
    end
end
