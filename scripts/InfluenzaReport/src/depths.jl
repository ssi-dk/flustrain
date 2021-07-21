struct Depths
    # depths: Depth where majority is not deletion
    depths::Vector{UInt32}
    # Cov: Where reference is not deletion
    refdepths::Vector{UInt32}
    mean_depth::Float64
    coverage::Float64
end

function Depths(depths::Vector{UInt32}, refdepths::Vector{UInt32})
    mean_depth = if isempty(depths)
        0.0
    else
        sum(depths) / length(depths)
    end
    coverage = if isempty(refdepths)
        0.0
    else
        sum(!iszero, refdepths) / length(refdepths)
    end
    Depths(depths, refdepths, mean_depth, coverage)
end

function load_depths(
    alnasms::Vector{SegmentTuple{Option{AlignedAssembly}}},
    depthspaths::Vector{String}
)::Vector{SegmentTuple{Option{Depths}}}

    res = SegmentTuple{Option{Depths}}[]
    # TODO: Parallelize?
    for (alnasm, depthspath) in zip(alnasms, depthspaths)
        mat = open(GzipDecompressorStream, depthspath) do io
            KMATools.parse_mat(io, depthspath)
        end
        depths = parse_depths(mat, alnasm, depthspath)
        add_depths_errors!(alnasm, depths)
        push!(res, depths)
    end
    return res
end

function parse_depths(
    mat::Vector{Tuple{String, Vector{Tuple{DNA, NTuple{6, UInt32}}}}},
    alnasm::SegmentTuple{Option{AlignedAssembly}},
    path::String
)::SegmentTuple{Option{Depths}}
    result = fill(none(Depths), N_SEGMENTS)
    for (header, rows) in mat
        segment::Segment = let
            p = findlast('_', header)
            s = p === nothing ? nothing : tryparse(Segment, strip(header[p+1:end]))
            s !== nothing ? s : error("Could not parse segment from file $path, header $(header)")
        end
        index = Integer(segment) + 0x01

        # If the alnasm is none, just leave depths as none, too
        if is_error(alnasm[index])
            continue
        end
        if !is_error(result[index])
            error("Segment $segment present twice in $string")
        end
        
        depthvec = UInt32[]
        refdepthvec = UInt32[]
        for (refbase, rowdepths) in rows
            rowdepth = sum(rowdepths[1:4]) # only add A,C,G,T, skip N
            if argmax(rowdepths) != 6
                push!(depthvec, rowdepth)
            end
            if !isgap(refbase)
                push!(refdepthvec, rowdepth)
            end
        end
        result[index] = some(Depths(depthvec, refdepthvec))
    end
    SegmentTuple(result)
end

function add_depths_errors!(
    alnasms::SegmentTuple{Option{AlignedAssembly}},
    depths::SegmentTuple{Option{Depths}}
)::Nothing
    for (m_alnasm, m_depth) in zip(alnasms, depths)
        alnasm = @unwrap_or m_alnasm continue
        depth = @unwrap_or m_depth continue

        if length(depth.depths) < 2*TERMINAL + 1
            push!(alnasm.errors, Influenza.ErrorTooShort(length(depth.depths)))
        end

        if depth.coverage < 0.9
            push!(alnasm.errors, Influenza.ErrorLowCoverage(depth.coverage))
        end

        n_lowdepth = count(i -> i < 25, @view depth.depths[TERMINAL + 1: end - TERMINAL])
        if !iszero(n_lowdepth)
            push!(alnasm.errors, Influenza.ErrorLowDepthBases(n_lowdepth))
        end
    end
end

function plot_depths(
    dir::AbstractString,
    basenames::Vector{String},
    depthsvec::Vector{SegmentTuple{Option{Depths}}}
)
    isdir(dir) || mkdir(dir)
    for (basename, depths) in zip(basenames, depthsvec)
        plt = plot(ylabel="Log10 depths", xticks=nothing, ylim=(-0.1, 5))
        for (index, m_depth) in enumerate(depths)
            depth = @unwrap_or m_depth continue
            segment = Segment(index - 1)
            ys = log10.(depth.depths)
            xs = range(0.0, stop=1.0, length=length(ys))
            plot!(plt, xs, ys, label=string(segment), legend=:outertopright)
        end
        savefig(plt, joinpath(dir, basename * ".pdf"))
    end
end