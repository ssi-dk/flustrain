using Plots
using CodecZlib

function get_blocks(path::AbstractString)::Vector{Vector{<:AbstractString}}
    io = open(path)
    lines = map(strip, eachline(GzipDecompressorStream(io)))
    close(io)
    empties = findall(isempty, lines)
    
    blocks = Vector{Vector{String}}()
    start = 1
    for i in empties
        push!(blocks, lines[start:i-1])
        start = i + 1
    end
    push!(blocks[start:end-1])
    blocks
end

function get_depths(path::AbstractString)::Vector{Tuple{String, Vector{UInt32}}}
    # Get a list of "blocks", each representing a segment
    blocks = get_blocks(path)
    result = Vector{Tuple{String, Vector{UInt32}}}()
    for block in blocks
        depths = UInt32[]
        segment = String(first(block)[2:end]) # looks like "#HA"
        for fields in (split(line) for line in @view(block[2:end]))
            first(first(fields)) == '-' && continue

            # A, C, G, T and also -, but not N
            depth = sum(x -> parse(UInt32, x), @view fields[2:5]) + parse(UInt32, fields[7])
            push!(depths, depth)
        end
        push!(result, (segment, depths))
    end
    result
end

function plot_depths!(plt::Plots.Plot, label::AbstractString, depth::Vector{UInt32})
    ys = log10.(depth)
    xs = range(0.0, stop=1.0, length=length(ys))
    plot!(plt, xs, ys, label=label, legend=:outertopright)
end

function main(outpath::AbstractString, inpath::AbstractString)
    plt = plot(xlabel="Nucleotide position", ylabel="Log10 depth")
    for (label, depth) in get_depths(inpath)
        plot_depths!(plt, label, depth)
    end

    savefig(plt, outpath)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], ARGS[2])
end
