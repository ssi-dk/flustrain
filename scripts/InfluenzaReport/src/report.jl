# Combine assembly and reference to a single object
struct ReferenceAssembly
    segment::Segment
    refseq::LongDNASeq
    asmseq::LongDNASeq
    insignificant::BitVector
    aln::PairwiseAlignment{LongDNASeq, LongDNASeq}
    aln_id::Float64
    orfdata::Vector{ORFData}
    messages::Vector{ErrorMessage}
    passed::Base.RefValue{Option{Bool}} # this is set at quality control
end

function ReferenceAssembly(ref::Reference, asm::Assembly)
    @assert ref.segment == asm.segment

    # For optimization: The large majority of time is spent on this alignment
    aln = pairalign(OverlapAlignment(), asm.seq, ref.seq, DEFAULT_DNA_ALN_MODEL).aln
    identity = unwrap(alignment_identity(aln))
    orfdata = map(ref.proteins) do protein
        ORFData(protein, aln, ref)
    end

    # Check for various errors here
    messages = ErrorMessage[]
    
    # Low identity to reference
    identity < 0.9 && push!(messages, ErrorMessage(important, "Identity to reference less than 90%"))

    # Insignificant bases
    n_insignificant = count(asm.insignificant)
    if !iszero(n_insignificant)
        severity = n_insignificant > 4 ? important : trivial
        push!(messages, ErrorMessage(severity, "$n_insignificant bases insignificantly called"))
    end

    # Ambiguous bases
    n_amb = count(isambiguous, asm.seq)
    if !iszero(n_amb)
        severity = n_amb > 4 ? important : trivial
        push!(messages, ErrorMessage(severity, "$n_amb central bases insignificantly called"))
    end

    ReferenceAssembly(asm.segment, ref.seq, asm.seq, asm.insignificant, aln,
        identity, orfdata, messages, Ref(none(Bool)))
end

"This creates a FASTA record with lowercase letters on uncertain positions"
function asm_dna_record(x::ReferenceAssembly, header::String)
    record = FASTA.Record(header, x.asmseq)
    data = record.data
    @assert length(record.sequence) == length(x.insignificant)
    for (isbad, seqpos) in zip(x.insignificant, record.sequence)
        if isbad
            data[seqpos] += 0x20 # this sets it to lowercase for ASCII
        end
    end
    return record
end

function segment_report_lines(segment::Segment, maybe_refasm::Option{ReferenceAssembly},
    extra::AbstractDict{String, Any}
)::Vector{String}
    beginning = rpad(string(segment) * ":", 4)

    # Check for presence
    if is_error(maybe_refasm)
        return ["FAIL " * beginning * " ERROR: Missing segment"]
    end
    refasm = unwrap(maybe_refasm)

    # Check segment length and return early if way too short
    if length(refasm.asmseq) < 2 * TERMINAL + 1
        refasm.passed[] = some(false)
        return [beginning * " ERROR: Sequence or depths length: $(length(efasm.asmseq))"]
    end
    
    # Add identity to the header line
    is_ok = true
    lines = String[beginning * " identity $(@sprintf "%.3f" (refasm.aln_id))"]

    # Depths contains Option{Tuple{String, Vector{ErrorMessage}}}, where the initial
    # string should be added to the header line.
    # We also put it up top, just because it's here anyway
    if haskey(extra, "depths")
        (depths_str, errors) = error_report(extra["depths"])
        lines[1] *= (' ' * depths_str)
        for message in errors
            push!(lines, format(message))
            is_ok &= is_trivial(message)
        end
    end

    # Add own messages from the ReferenceAssembly itself
    for message in refasm.messages
        push!(lines, format(message))
        is_ok &= is_trivial(message)
    end

    # Add the errors of the ORFdata
    for orfdata in refasm.orfdata, message in orfdata.errors
        push!(lines, format(message))
        is_ok &= is_trivial(message)
    end

    # Now add all the rest (extra checks and data)
    for (name, messages) in extra
        name == "depths" && continue
        for message in messages
            push!(lines, format(message))
            is_ok &= is_trivial(message)
        end
    end

    # Finally, we add the OK status to header of lines
    lines[1] = (is_ok ? "     " : "FAIL ") * lines[1]
    refasm.passed[] = some(is_ok)
    return lines
end

function basename_report_lines(
    maybe_refasm_tuple::SegmentTuple{Option{ReferenceAssembly}},
    extra_tuple::SegmentTuple{<:AbstractDict},
)::SegmentTuple{Vector{String}}
    result = Vector{Vector{String}}()
    for (segment, maybe_refasm, extra) in zip(instances(Segment), maybe_refasm_tuple, extra_tuple)
        lines = segment_report_lines(segment, maybe_refasm, extra)

        # Indent all lines except first
        for i in 2:lastindex(lines)
            lines[i] = '\t' * lines[i]
        end

        push!(result, lines)
    end
    return SegmentTuple(result)
end

function write_report(io::IO, basenames::Vector{<:AbstractString},
    basename_reports::Vector{SegmentTuple{Vector{String}}}
)
    for (basename, reports) in zip(basenames, basename_reports)
        println(io, basename)
        for lines in reports
            for line in lines
                println(io, '\t', line)
            end
        end
        print(io, '\n')
    end
end

function write_consensus(dirname::String, basenames::Vector{String},
    maybe_refasm_tuples::Vector{SegmentTuple{Option{ReferenceAssembly}}}
)
    isdir(dirname) || mkdir(dirname)
    for (basename, maybe_refasm_tuple) in zip(basenames, maybe_refasm_tuples)
        subdir = joinpath(dirname, basename)
        isdir(subdir) || mkdir(subdir)
        cons_dna_writer = open(FASTA.Writer, joinpath(subdir, "consensus.fna"))
        cons_aa_writer = open(FASTA.Writer, joinpath(subdir, "consensus.faa"))
        cura_dna_writer = open(FASTA.Writer, joinpath(subdir, "curated.fna"))
        cura_aa_writer = open(FASTA.Writer, joinpath(subdir, "curated.faa"))
        for maybe_refasm in maybe_refasm_tuple
            refasm = @unwrap_or maybe_refasm continue
            is_ok = unwrap(refasm.passed[])
            dna_record = asm_dna_record(refasm, basename * '_' * string(refasm.segment))
            write(cons_dna_writer, dna_record)
            is_ok && write(cura_dna_writer, dna_record)
            for orfdata in refasm.orfdata
                header = basename * '_' * string(orfdata.variant) * '_'
                header *= join(["$(first(i))-$(last(i))" for i in orfdata.orfs], ',')
                record = FASTA.Record(header, orfdata.aaseq)
                write(cons_aa_writer, record)
                is_ok && write(cura_aa_writer, record)
            end
        end
        close(cons_dna_writer)
        close(cons_aa_writer)
        close(cura_dna_writer)
        close(cura_aa_writer)
    end
end

# This function's API is meant for snakemake
function illumina_snakemake_entrypoint(
    report_path::AbstractString, # output report
    ref_dir::AbstractString, # dir of .fna + .jls ref files
    aln_dir::AbstractString, # dir of medaka / kma aln
    cons_dir::AbstractString,
    depths_plot_dir::AbstractString
)::Nothing
    basenames = sort!(readdir(aln_dir))

    # Load assemblies
    maybe_asm_tuples = basenames |> Map() do basename
        load_assembly(joinpath(aln_dir, basename, "kma2.fsa"), true)
    end |> Folds.collect

    # Get references
    maybe_ref_tuples = load_references(maybe_asm_tuples, ref_dir)

    # Merge assemblies and references to maybe_refasms.
    maybe_refasm_tuples = zip(maybe_asm_tuples, maybe_ref_tuples) |> 
    Map() do (maybe_asm_tuple, maybe_ref_tuple)
        ntuple = SegmentTuple(zip(maybe_asm_tuple, maybe_ref_tuple))
        map(ntuple) do (maybe_asm, maybe_ref)
            (is_error(maybe_asm) || is_error(maybe_ref)) && return none(ReferenceAssembly)
            some(ReferenceAssembly(unwrap(maybe_ref), unwrap(maybe_asm)))
        end
    end |> Folds.collect

    # Extra data:
    extras = [ntuple(i -> Dict{String, Any}(), N_SEGMENTS) for i in basenames]

    # Extra: Add depth (and plot it!)
    isdir(depths_plot_dir) || mkdir(depths_plot_dir)
    for (basename, dict_tuple) in zip(basenames, extras)
        depths_matrices = from_mat(joinpath(aln_dir, basename, "kma1.mat.gz"))
        depths = map(depths_matrices) do maybe_matrix
            and_then(depths_from_mat, Depths, maybe_matrix)
        end

        # Plot depths
        plot_depths(joinpath(depths_plot_dir, basename * ".pdf"), depths)

        # Add depth info to extra dict
        for (dict, maybe_depth) in zip(dict_tuple, depths)
            if !is_error(maybe_depth)
                dict["depths"] = unwrap(maybe_depth)
            end
        end

        # Add in information for superinfection
        for (dict, maybe_matrix) in zip(dict_tuple, depths_matrices)
            if !is_error(maybe_matrix)
                if is_superinfected(unwrap(maybe_matrix))
                    dict["superinfection"] = [ErrorMessage(trivial, "[ POSSIBLE SUPERINFECTION ]")]
                end
            end
        end
    end

    # Extra: Add kma2 identity check
    for (basename, dict_tuple) in zip(basenames, extras)
        second_kma_file = joinpath(aln_dir, basename, "kma2.res")
        for (dict, maybe_identity_error) in zip(dict_tuple, check_second_kma_file(second_kma_file))
            if !is_error(maybe_identity_error)
                dict["second_identity"] = [unwrap(maybe_identity_error)]
            end
        end
    end

    # Get segment report lines
    basename_reports = zip(maybe_refasm_tuples, extras) |> Map() do (maybe_refasm_tup, extra_tup)
        basename_report_lines(maybe_refasm_tup, extra_tup)
    end |> Folds.collect

    # Create report itself
    open(report_path, "w") do io
        write_report(io, basenames, basename_reports)
    end
    
    # Write out consensus
    write_consensus(cons_dir, basenames, maybe_refasm_tuples)
    return nothing
end

function nanopore_snakemake_entrypoint(
    report_path::AbstractString, # output report
    ref_dir::AbstractString, # dir of .fna + .jls ref files
    aln_dir::AbstractString, # dir of medaka / kma aln
    cons_dir::AbstractString,
)::Nothing
    basenames = sort!(readdir(aln_dir))

    # Load assemblies
    maybe_asm_tuples = basenames |> Map() do basename
        load_assembly(joinpath(aln_dir, basename, "medaka/consensus.fasta"), false)
    end |> Folds.collect

    # Get references
    maybe_ref_tuples = load_references(maybe_asm_tuples, ref_dir)

    # Merge assemblies and references to maybe_refasms.
    maybe_refasm_tuples = zip(maybe_asm_tuples, maybe_ref_tuples) |> 
    Map() do (maybe_asm_tuple, maybe_ref_tuple)
        ntuple = SegmentTuple(zip(maybe_asm_tuple, maybe_ref_tuple))
        map(ntuple) do (maybe_asm, maybe_ref)
            (is_error(maybe_asm) || is_error(maybe_ref)) && return none(ReferenceAssembly)
            some(ReferenceAssembly(unwrap(maybe_ref), unwrap(maybe_asm)))
        end
    end |> Folds.collect

    # Extra data:
    extras = [ntuple(i -> Dict{String, Any}(), N_SEGMENTS) for i in basenames]

    # Get segment report lines
    basename_reports = zip(maybe_refasm_tuples, extras) |> Map() do (maybe_refasm_tup, extra_tup)
        basename_report_lines(maybe_refasm_tup, extra_tup)
    end |> Folds.collect

    # Create report itself
    open(report_path, "w") do io
        write_report(io, basenames, basename_reports)
    end
    
    # Write out consensus
    write_consensus(cons_dir, basenames, maybe_refasm_tuples)
    return nothing
end