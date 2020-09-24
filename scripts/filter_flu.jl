using CodecZlib
using FASTX
using SHA
using Comonicon

const ENDINGS = ("fastq.gz", ".fq.gz")

struct Hash # 128 bits is enough
    a::UInt64
    b::UInt64
end

function Hash(st::Union{AbstractString, AbstractArray{UInt8}})
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
"Parse the out.frag.gz file and extract a set of mapping readnames"
function build_present_set(path::AbstractString, minscore::Int)::Set{Hash}
    open(path) do file
        io = GzipDecompressorStream(file)
        s = sizehint!(Set{Hash}(), 100000)
        for (lineno, line) in enumerate(eachline(io))
            fields = split(line)
            if length(fields) != 7
                throw(ArgumentError("Line $lineno of file $path does not contain 7 fields"))
            end
            score = parse(Int, fields[3])
            score > minscore && push!(s, Hash(last(split(line))))
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

@main function filter_flu(dstdir, srcdir, kmadb, minscore::Int=100)
    for i in (dstdir, srcdir)
        isdir(i) || throw(SystemError("Directory not found: $i"))
    end
        
    mkdir(dstdir)
    pairs = get_fastqs(srcdir)
    Threads.@threads for (fw, rv) in pairs
        filter_fastq(dstdir, srcdir, fw, rv, db)
    end
end
