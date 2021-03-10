# Here, we load template identities from kma2.res. If a template has <99.5% identity,
# it means the first assembly did not converge.
function load_kma_file(resfilename::AbstractString)::SegmentTuple{Option{Float64}}
    open(resfilename) do io
        fields = Vector{SubString{String}}(undef, 11)
        lines = eachline(io)
        header, _ = iterate(lines)::NTuple{2, Any}
        result = fill(none(Float64), length(instances(Segment)))
        @assert header == "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value"
        for fields in lines |> Map(strip) ⨟ Filter(!isempty) ⨟ Map(x -> split!(fields, line, UInt8('\t')))
            segment = unwrap(parse(Segment, strip(last(rsplit(first(fields), '_', limit=2)))))
            segment_index = reinterpret(UInt8, segment) + 0x01
            is_error(result[segment_index]) || error("Segment $(string(segment)) present twice in file $resfilename")
            result[segment_index] = some(parse(Float64, fields[5]) / 100)
        end
        SegmentTuple(result)
    end
end

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
end

function ReferenceAssembly(ref::Reference, asm::Assembly)
    @assert ref.segment == asm.segment
    aln = pairalign(OverlapAlignment(), asm.seq, ref.seq, ALN_MODEL).aln
    identity = unwrap(get_identity(aln))
    orfdata = map(ref.proteins) do protein
        ORFData(protein, aln, ref)
    end

    # Check for various errors here
    messages = ErrorMessage[]
    
    # Low identity to reference
    identity < 0.9 && push!(messages, ErrorMessage(important, "Identity less than 90%"))

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
        identity, orfdata, messages)
end

# To do: Add other error messages from e.g. depths
function segment_report_lines(segment::Segment, maybe_refasm::Option{ReferenceAssembly},
    extra::NamedTuple{S, <:Tuple{Vararg{<:Option}}} where S
)::Tuple{Bool, Vector{String}}
    beginning = rpad(string(segment) * ":", 4)

    # Check for presence
    if is_error(maybe_refasm)
        return (false, [beginning * " ERROR: Missing segment"])
    end
    refasm = unwrap(maybe_refasm)

    # Check segment length and return early if way too short
    if length(refasm.asmseq) < 2 * TERMINAL + 1
        return (false, [beginning * " ERROR: Sequence or depths length: $(length(efasm.asmseq))"])
    end
    
    # Add identity to the header line
    is_ok = true
    lines = String[beginning * " identity $(@sprintf "%.3f" (refasm.aln_id))"]

    # Depths contains Option{Tuple{String, Vector{ErrorMessage}}}, where the initial
    # string should be added to the header line.
    # We also put it up top, just because it's here anyway
    if haskey(extra, :depths) && !is_error(extra.depths)
        (depths_str::String, errs::Vector{ErrorMessage}) = unwrap(extra.depths)
        lines[1] *= (' ' * depths_str)
        for message in errs
            push!(lines, format(message))
            is_ok &= is_trivial(message)
        end
    end

    # Add own messages
    for message in refasm.messages
        push!(lines, format(message))
        is_ok &= is_trivial(message)
    end

    # Add the errors of the ORFdata
    for orfdata in refasm.orfdata, message in orfdata.errors
        push!(lines, format(message))
        is_ok &= is_trivial(message)
    end

    # Now add all the rest.
    for (symbol, maybe_msgs) in pairs(extra)
        symbol == :depths && continue
        msgs = @unwrap_or maybe_msgs continue
        for msg in msgs
            push!(lines, format(msg))
            is_ok &= is_trivial(message)
        end
    end
 
    return (is_ok, lines)
end

function basename_report_lines(
    maybe_refasm_tuple::SegmentTuple{Option{ReferenceAssembly}},
    extra_tuple::SegmentTuple{NamedTuple},
)::Tuple{SegmentTuple{Bool}, SegmentTuple{Vector{String}}}
    is_oks = Bool[]
    line_vecs = Vector{String}[]
    for (maybe_refasm, extra) in zip(maybe_refasm_tuple, extra_tuple)
        is_ok, lines = segment_report_lines(maybe_refasm, extra)

        # Indent all lines except first
        for i in 2:lastindex(lines)
            lines[i] = '\t' * lines[i]
        end

        push!(is_oks, is_ok)
        push!(line_vecs, lines)
    end
    (SegmentTuple(is_oks), SegmentTuple(line_vecs))
end

# This function's API is meant for snakemake
function illumina_snakemake_entrypoint(
    report_path::AbstractString, # output report
    ref_dir::AbstractString, # dir of .fna + .jls ref files
    aln_dir::AbstractString, # dir of medaka / kma aln
    depths_plot_dir::AbstractString
)
    basenames = sort!(readdir(aln_dir))

    # Load assemblies
    maybe_asm_tuples = basenames |> Map() do basename
        load_assembly(joinpath(aln_dir, basename, "kma2.fsa"), true)
    end |> Folds.collect

    # Get references
    maybe_ref_tuples = load_references(maybe_asm_tuples, ref_dir)

    # Convert to maybe_refasms.
    maybe_refasm_tuples = zip(maybe_asm_tuples, maybe_ref_tuples) |> 
    Map() do (maybe_asm_tuple, maybe_ref_tuple)
        ntuple = SegmentTuple(zip(maybe_asm_tuple, maybe_ref_tuple))
        map(ntuple) do (maybe_asm, maybe_ref)
            (is_error(maybe_asm) || is_error(maybe_ref)) && return none(ReferenceAssembly)
            some(ReferenceAssembly(unwrap(maybe_ref), unwrap(maybe_asm)))
        end
    end |> Folds.collect

    # Load depths
    maybe_depths_tuples = basenames |> Map() do basename
        from_mat(joinpath(aln_dir, basename, "kma1.mat.gz"))
    end |> Folds.collect

    # Extras
    #extras_tuples = zip(maybe_depths_tuples,) |> Map() do (m_depth)


    (maybe_refasm_tuples, maybe_depths_tuples)
end
