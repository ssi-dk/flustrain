function illumina_snakemake_entrypoint(
    report_path::AbstractString, # output report
    ref_dir::AbstractString, # dir of .fna + .jls ref files
    aln_dir::AbstractString, # dir of medaka / kma aln
    cons_dir::AbstractString,
    depths_plot_dir::AbstractString
)::Nothing
    basenames = sort!(readdir(aln_dir))

    asm_paths = [joinpath(aln_dir, basename, "kma2.fsa") for basename in basenames]
    aln_asms = load_aligned_assemblies(asm_paths, ref_dir, true)

    kma2_paths = [joinpath(aln_dir, basename, "kma2.res") for basename in basenames]
    kma2_identity_check(aln_asms, kma2_paths)

    depthpaths = [joinpath(aln_dir, basename, "kma1.mat.gz") for basename in basenames]
    depths = load_depths(aln_asms, depthpaths)
    plot_depths(depths_plot_dir, basenames, depths)

    passes = report(report_path, basenames, aln_asms, depths)
    write_consensus(cons_dir, basenames, aln_asms, passes)
end

function nanopore_snakemake_entrypoint(
    report_path::AbstractString, # output report
    ref_dir::AbstractString, # dir of .fna + .jls ref files
    aln_dir::AbstractString, # dir of medaka / kma aln
    cons_dir::AbstractString,
    depths_plot_dir::AbstractString
)::Nothing
    basenames = sort!(readdir(aln_dir))

    asm_paths = [joinpath(aln_dir, basename, "medaka", "consensus.fasta") for basename in basenames]
    aln_asms = load_aligned_assemblies(asm_paths, ref_dir, false)

    depthpaths = [joinpath(aln_dir, basename, "kma1.mat.gz") for basename in basenames]
    depths = load_depths(aln_asms, depthpaths)
    plot_depths(depths_plot_dir, basenames, depths)

    passes = report(report_path, basenames, aln_asms, depths)
    write_consensus(cons_dir, basenames, aln_asms, passes)
end

# This function has a bit weird logical flow, because it's easiest to print
# at the same time as determining whether the segment has passed.
# But we also need to prefix failed segments with "FAIL".
# So what we do is to print to an IOBuffer, then prefix any segment with FAIL.
function report(
    report_path::AbstractString,
    basenames::Vector{<:AbstractString},
    alnasms::Vector{SegmentTuple{Option{AlignedAssembly}}},
    depths::Vector{SegmentTuple{Option{Depths}}},
)::Vector{SegmentTuple{Bool}}
    result = SegmentTuple{Bool}[]
    open(report_path, "w") do io
        for (basename, alnasm, depth) in zip(basenames, alnasms, depths)
            push!(result, report(io, basename, alnasm, depth))
        end
        print(io, '\n')
    end
    return result
end

function report(
    io::IO,
    basename::AbstractString,
    alnasms::SegmentTuple{Option{AlignedAssembly}},
    depths::SegmentTuple{Option{Depths}}
)::SegmentTuple{Bool}
    passes = Bool[]
    println(io, basename)
    for (i, (alnasm, depth)) in enumerate(zip(alnasms, depths))
        segment = Segment(i - 1)
        (buf, pass) = report(alnasm, depth)
        print(io,
            '\t', pass ? "     " : "FAIL ",
            rpad(string(segment) * ":", 4),
        )
        write(io, take!(buf))
        push!(passes, pass)
    end
    print(io, '\n')
    return SegmentTuple(passes)
end

function report(
    maybe_alnasm::Option{AlignedAssembly},
    maybe_depth::Option{Depths},
)::Tuple{IOBuffer, Bool}
    buf = IOBuffer()
    passed = true
    if is_error(maybe_alnasm)
        println(buf, " Missing segment")
        return (buf, false)
    end
    (alnasm, depth) = (unwrap(maybe_alnasm), unwrap(maybe_depth))
    print_segment_header(buf, alnasm, depth)
    for error in alnasm.errors
        p = pass(error)
        println(buf, "\t\t", (p ? "      " : "ERROR "), error)
        passed &= p
    end

    IndelError = Union{Influenza.ErrorFrameShift, Influenza.ErrorIndelTooBig}
    for protein in alnasm.proteins

        # There can be a ton of indel errors. We're not interested in printing all of them,
        # so if there are a lot of indel errors, we just replace it with a ErrorTooManyIndels.
        if count(i -> i isa IndelError, protein.errors) > 3
            N = length(protein.errors)
            filter!(i -> !isa(i, IndelError), protein.errors)
            n_indel_errors = length(protein.errors) - N
            push!(protein.errors, ErrorTooManyIndels(n_indel_errors))
        end

        for error in protein.errors
            p = pass(error)
            println(buf, "\t\t", (p ? "      " : "ERROR "), protein.variant, ": ", error)
            passed &= p
        end
    end
    return (buf, passed)
end

function print_segment_header(buf::IOBuffer, alnasm::AlignedAssembly, depth::Depths)
    print(buf, " identity $(@sprintf "%.3f" (alnasm.identity))")
    println(buf, " depth $(@sprintf "%.2e" depth.mean_depth) coverage $(@sprintf "%.3f" depth.coverage)")
end

# fallback: segment errors fails, proteinerrors fail if protein is critical
pass(::Influenza.InfluenzaError) = false
pass(::Influenza.ProteinError, protein::Protein) = !is_important(protein)
pass(x::Influenza.ErrorInsignificant) = x.n_insignificant < 5
pass(x::Influenza.ErrorAmbiguous) = x.n_ambiguous < 5

# Fail if 15 or more AAs missing
function pass(x::Influenza.ErrorEarlyStop, protein::Protein)
    !(is_important(protein) && x.expected_naa - x.observed_naa > 14)
end

# Not critical
pass(x::Influenza.ErrorLateStop, protein::Protein) = true
