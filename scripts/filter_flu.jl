using CodecZlib
using FASTX
using SHA
using Comonicon

const ENDINGS = ("fastq.gz", ".fq.gz")

struct Hash # 128 bits is enough
    a::UInt64
    b::UInt64
end

function Hash(rec::FASTQ.Record)
    isempty(rec.identifier) && throw(ArgumentError("Read empty"))
    sh = sha256(view(rec.data, rec.identifier))
    GC.@preserve sh begin
        p = Ptr{UInt64}(pointer(sh))
        h = Hash(unsafe_load(p, 1), unsafe_load(p, 2))
    end
    return h
end

function Hash(st::AbstractString)
    sh = sha256(st)
    GC.@preserve st begin
        p = Ptr{UInt64}(pointer(sh))    
        h = Hash(unsafe_load(p, 1), unsafe_load(p, 2))
    end
    return h
end

##
function parse_filename(s::AbstractString)
    m = match(r"^([^_]+_S\d+)_L001_R(1|2)_001.fastq.gz$", s)
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

##
function build_present_set(path::AbstractString)::Set{Hash}
    open(path) do file
        io = GzipDecompressorStream(file)
        s = sizehint!(Set{Hash}(), 100000)
        for line in eachline(io)
            push!(s, Hash(last(split(line))))
        end
        s
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
    run(`kma -ipe $fw $rv -o $dir/kmaout -t_db $db`)
end

function filter_fastq(dstdir, srcdir, fw, rv, db)
    fwin = joinpath(srcdir, fw)
    fwout = joinpath(dstdir, fw)
    rvin = joinpath(srcdir, rv)
    rvout = joinpath(dstdir, rv)

    hashset = mktempdir(".") do dir
        kma_map(dir, fwin, rvin, db)
        build_present_set(joinpath(dir, "kmaout.frag.gz"))
    end

    filter_fastq(hashset, fwin, fwout)
    filter_fastq(hashset, rvin, rvout)
end

@main function filter_flu(dstdir, srcdir, db)
    mkdir(dstdir)
    pairs = get_fastqs(srcdir)
    Threads.@threads for (fw, rv) in pairs
        filter_fastq(dstdir, srcdir, fw, rv, db)
    end
end
