struct Depths
    depths::Vector{UInt32}
end

function error_report(maybe_depths::Option{Depths})::Option{Tuple{String, Vector{ErrorMessage}}}
    depths = @? maybe_depths
    if length(depths.depths) < 2*TERMINAL+1
        return ("ERROR: DEPTHS TOO SHORT", ErrorMessage[])
    end
    mean_depth = sum(UInt, depths.depths) / length(depths.depths)
    coverage = count(!iszero, depths.depths) / length(depths.depths)
    header = "depth $(@sprintf "%.2e" mean_depth) coverage $(@sprintf "%.3f" coverage)"
    
    messages = ErrorMessage[]
    
    # Low coverage
    if coverage < 0.9 &&
        push!(messages, ErrorMessage(important, "Coverage is less than 90%"))
    end

    # Low depth
    lowdepth = count(x -> x < 25, @view depths.depths[TERMINAL + 1: end - TERMINAL])
    if !iszero(lowdepth)
        severity = lowdepth > 4 ? important : trivial
        push!(messages, ErrorMessage(severity, "$lowdepth central bases with depth < 25"))
    end

    some((header, messages))
end

"Given a path to a .mat.gz file, return a dict of an optional vec of depths"
function from_mat(matpath::AbstractString)::SegmentTuple{Option{Depths}}
    open(matpath) do io
        result = fill(none(Vector{UInt32}), length(instances(Segment)))
        segment = nothing
        depths = UInt32[]
        fields = Vector{SubString{String}}(undef, 7)
        linedepths = Vector{UInt32}(undef, 6)
        for line in eachline(GzipDecompressorStream(io))
            if isempty(line)
                if !isempty(depths)
                    segment_index = reinterpret(UInt8, segment::Segment) + 0x01
                    is_error(result[segment_index]) || error("Segment $segment present twice in $matpath")
                    result[segment_index] = some(depths)
                    depths = UInt32[]
                end
                continue
            end
            if startswith(line, '#')
                # Headers look like "HEADER_HA"
                p = findlast(isequal('_'), line)
                p === nothing && error("Found header \"$(line)\", expected pattern \"HEADER_SEGMENT\"")
                segment = unwrap(parse(Segment, line[p+1:end]))
                continue
            end
            split!(fields, line, UInt8('\t'))
            @inbounds for i in 1:6
                linedepths[i] = parse(UInt32, fields[i+1])
            end
            # Skip lines with majority as deletion
            argmax(linedepths) == 6 && continue
            total_depth = sum(linedepths)
            # We don't include N in depths!
            depth = total_depth - linedepths[5]
            push!(depths, depth)
        end
        SegmentTuple(map(result) do maybe_vector
            and_then(Depths, Depths, maybe_vector)
        end)
    end
end

function from_samtools_depth(path::AbstractString)::SegmentTuple{Option{Depths}}
    open(path) do io
        position = 1
        result = fill(none(Vector{UInt32}), length(instances(Segment)))
        fields = Vector{SubString{String}}(undef, 3)
        for line in eachline(GzipDecompressorStream(io))
            split!(fields, line, UInt8('\t'))
            segment = unwrap(parse(Segment, first(fields)[findlast('_', first(fields)) + 1:end]))
            segment_index = reinterpret(UInt8, segment) + 0x01
            vector = if is_error(result[segment_index])
                position = 1
                v = UInt32[]
                result[segment_index] = some(v)
                v
            else
                position += 1
                unwrap(result[segment_index])
            end
            @assert length(vector) == position - 1
            @assert position == parse(UInt, fields[2])
            push!(vector, parse(UInt32, fields[3]))
        end
        SegmentTuple(map(result) do maybe_vector
            and_then(Depths, Depths, maybe_vector)
        end)
    end
end

# This function is thread-unsafe, and must be called in serial.
function plot_depths(path::AbstractString, depths::SegmentTuple{Option{Depths}})
    plt = plot(ylabel="Log10 depths", xticks=nothing, ylim=(-0.1, 5))
    for (index, maybe_depths) in depths
        data = @unwrap_or maybe_depths continue
        segment = Segment(index - 1)
        ys = log10.(data.depths)
        xs = range(0.0, stop=1.0, length=length(ys))
        plot!(plt, xs, ys, label=string(segment), legend=:outertopright)
    end
    savefig(plt, path)
end