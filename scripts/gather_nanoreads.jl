using CodecZlib

const SKIPPABLE = ["unclassified"]

function copy_compress(outpath::String, indir::String)
    isdir(indir) || error("Not a directory: $indir")
    buffer = Vector{UInt8}(undef, 2^16)
    io_out = GzipCompressorStream(open(outpath, "w"))
    for file in readdir(indir)
        infile = open(joinpath(indir, file))
        io_in = endswith(file, ".gz") ? GzipDecompressorStream(infile) : infile
        while true
            nbytes = readbytes!(io_in, buffer)
            iszero(nbytes) && break
            unsafe_write(io_out, pointer(buffer), nbytes)
        end
        close(io_in)
    end
    close(io_out)
end

function main(indir, outdir)
    isdir(indir) || error("Not a directory: $indir")
    mkdir(outdir)

    subdirs = setdiff(iterate(walkdir(indir))[1][2], SKIPPABLE)
    tasks = Task[]
    for subdir in subdirs
        inpath = joinpath(indir, subdir)
        outpath = joinpath(outdir, subdir * ".fastq.gz")
        push!(tasks, Base.Threads.@spawn copy_compress(outpath, inpath))
    end
    foreach(wait, tasks)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        println("Usage: julia gather_nanoreads.jl inpath outpath")
        exit(1)
    else
        main(ARGS...)
    end
end