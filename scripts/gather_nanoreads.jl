const SKIPPABLE = ["unclassified"]

function main(indir, outdir)
    isdir(indir) || error("Not a directory: $indir")
    mkdir(outdir)

    # Load files
    buffer = Vector{UInt8}(undef, 2^16)
    for dir in setdiff(iterate(walkdir(indir))[1][2], SKIPPABLE)
        open(joinpath(outdir, dir * ".fastq"), "w") do outfile
            for file in readdir(joinpath(indir, dir))
                open(joinpath(indir, dir, file)) do infile
                    while true
                        nbytes = readbytes!(infile, buffer)
                        iszero(nbytes) && break
                        unsafe_write(outfile, pointer(buffer), nbytes)
                    end
                end
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        println("Usage: julia gather_nanoreads.jl inpath outpath")
        exit(1)
    else
        main(ARGS...)
    end
end