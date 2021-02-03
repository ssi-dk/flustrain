using Plots
using CodecZlib

function get_depth(path::AbstractString)
    # If the file is empty (e.g. it was created by a touch op), 
    # just return a small list of zeros
    iszero(stat(path).size) && return zeros(UInt32, 10)
    
    depths = UInt32[]
    open(path) do io
        lines = eachline(GzipDecompressorStream(io))
        iterate(lines) # skip header
        stripped = (rstrip(line) for line in lines)
        for fields in (split(strip) for strip in stripped if !isempty(strip))
            first(first(fields)) == '-' && continue

            # A, C, G, T and also -, but not N
            depth = sum(x -> parse(UInt32, x), @view fields[2:5]) + parse(UInt32, fields[7])
            push!(depths, depth)
        end
        depths
    end
end

function plot_depths!(plt::Plots.Plot, inpath::AbstractString, label::AbstractString)
    ys = log10.(get_depth(inpath))
    xs = range(0.0, stop=1.0, length=length(ys))
    plot!(plt, xs, ys, label=label, legend=:outertopright)
end

function main(outpath::AbstractString, inpaths::Vector{String})
    segments = [first(split(basename(x), ".", limit=2)) for x in inpaths]
    plt = plot(xlabel="Nucleotide position", ylabel="Log10 depth")
    for (inpath, label) in zip(inpaths, segments)
        plot_depths!(plt, inpath, label)
    end

    savefig(plt, outpath)
end

main(ARGS[1], ARGS[2:end])
