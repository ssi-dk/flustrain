struct Depths
    depths::Vector{UInt32}
end

function error_report(depths::Depths)::Tuple{String, Vector{ErrorMessage}}
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

    return (header, messages)
end

depths_from_mat(mat::Matrix{<:Integer}) = Depths(vec(sum(mat, dims=1)))

function from_mat(matpath::AbstractString)::SegmentTuple{Option{Matrix{UInt32}}}
    open(matpath) do io
        result = fill(none(Matrix{UInt32}), N_SEGMENTS)
        segment = nothing
        depths = UInt32[]
        fields = Vector{SubString{String}}(undef, 7)
        linedepths = Vector{UInt32}(undef, 6)
        for line in eachline(GzipDecompressorStream(io))
            if isempty(line)
                if !isempty(depths)
                    segment_index = reinterpret(UInt8, segment::Segment) + 0x01
                    is_error(result[segment_index]) || error("Segment $segment present twice in $matpath")
                    result[segment_index] = some(reshape(depths, 4, :))
                    depths = UInt32[]
                end
                continue
            end
            if startswith(line, '#')
                # Headers look like "HEADER_HA"
                p = findlast(isequal('_'), line)
                p === nothing && error("Found header \"$(line)\", expected pattern \"HEADER_SEGMENT\"")
                segment = tryparse(Segment, line[p+1:end])::Segment
                continue
            end
            split!(fields, line, UInt8('\t'))
            @inbounds for i in 1:6
                linedepths[i] = parse(UInt32, fields[i+1])
            end
            # Skip lines with majority as deletion
            argmax(linedepths) == 6 && continue

            # Only append A, C, G, T counts, skip N and gap
            append!(depths, view(linedepths, 1:4))
        end
        SegmentTuple(result)
    end
end

"""
This splits the matrix into the most common (top) allele and the second-most common (bottom).
It uses heuristics to try to not get too much fluctuations in the bottom's depth.
"""
function split_matrix(m::AbstractMatrix{<:Integer})
    top, bottom = Vector{eltype(m)}(undef, size(m, 2)), Vector{eltype(m)}(undef, size(m, 2))
    bottomi, topi = nothing, zero(eltype(m))
    delta_max = 2.0
    v = zeros(eltype(m), size(m, 2))
    for col in 1:size(m, 2)
        topy, bottomy = sort!(copyto!(v, view(m, :, col)), rev=true)
        if bottomi === nothing || max(bottomy, bottomi) / min(bottomy, bottomi) < delta_max
            bottomi = bottomy
            topi = topy
        end
        bottom[col] = bottomi
        top[col] = topi
        bottomi = ifelse(bottomi < 10, nothing, bottomi)
    end
    top, bottom
end

#=
function split_matrix(m::AbstractMatrix{<:Integer})
    @assert size(m, 1) == 4 "Must have 4 rows for nucleotides A, C, G and T."

    # A site is a "minority variant" if it is the 2nd most commonly called nucleotide,
    # it has a depth > 10, is at least 0.1 of the major variant's depth, and at least
    # 0.5 of last minor variant base depth.
    top_deps = Vector{eltype(m)}(undef, size(m, 2))
    bottom_deps = similar(top_dep)
    top_seq, bottom_seq = DNA[], DNA[]
    last_bottom_dep = nothing
    bases_since_bottom = 0
    indexbuffer = zeros(Int, 4)
    BASES = (DNA_A, DNA_C, DNA_G, DNA_T)
    for col in 1:size(m, 2)
        top_index, bottom_index = sortperm!(indexbuffer, view(m, :, col), rev=true)
        top_dep, bottom_dep = m[top_index, col], m[bottom_index, col]
        top_base, bottom_base = BASES[top_index], BASES[bottom_index]

        if (
            bottom_dep ≥ 10 &&
            bottom_dep ≥ 0.1 * top_dep &&
            last_bottom_dep === nothing || bottom_dep ≥ 0.5 * last_bottom_dep
        )
            bases_since_bottom = 0
            last_bottom_dep = bottom_dep
        else
            bases_since_bottom += 1
        end





    top, bottom = Vector{eltype(m)}(undef, size(m, 2)), Vector{eltype(m)}(undef, size(m, 2))
    bottomi, topi = nothing, zero(eltype(m))
    delta_max = 2.0
    v = zeros(eltype(m), size(m, 2))
    for col in 1:size(m, 2)
        topy, bottomy = sort!(copyto!(v, view(m, :, col)), rev=true)
        if bottomi === nothing || max(bottomy, bottomi) / min(bottomy, bottomi) < delta_max
            bottomi = bottomy
            topi = topy
        end
        bottom[col] = bottomi
        top[col] = topi
        bottomi = ifelse(bottomi < 10, nothing, bottomi)
    end
    top, bottom
end
=#

"This is far from certain - it just gives a hint. Takes the matrix from from_mat"
function is_superinfected(m::AbstractMatrix{<:Integer})
    top, bottom = split_matrix(m)
    topmean = sum(top) / length(top)
    topmean > 250 && topmean < (10 * sum(bottom) / length(bottom))
end

function from_samtools_depth(path::AbstractString)::SegmentTuple{Option{Depths}}
    open(path) do io
        position = 1
        result = fill(none(Vector{UInt32}), N_SEGMENTS)
        fields = Vector{SubString{String}}(undef, 3)
        for line in eachline(GzipDecompressorStream(io))
            split!(fields, line, UInt8('\t'))
            segment::Segment = tryparse(Segment, first(fields)[findlast('_', first(fields)) + 1:end])
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
    for (index, maybe_depths) in enumerate(depths)
        data = @unwrap_or maybe_depths continue
        segment = Segment(index - 1)
        ys = log10.(data.depths)
        xs = range(0.0, stop=1.0, length=length(ys))
        plot!(plt, xs, ys, label=string(segment), legend=:outertopright)
    end
    savefig(plt, path)
end