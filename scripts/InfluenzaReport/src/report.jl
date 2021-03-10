const TERMINAL = 25

# Here, we load template identities from kma2.res. If a template has <99.5% identity,
# it means the first assembly did not converge.
function load_kma_file(resfilename::AbstractString)::SegmentTuple{Option{Float64}}
    open(resfilename) do io
        fields = Vector{SubString{String}}(undef, 11)
        lines = eachline(io)
        header, _ = iterate(lines)::NTuple{2, Any}
        result = fill(none(Float64), length(instances(Segment)))
        @assert header == "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value"
        for fields in (lines |> imap(strip) |> ifilter(!isempty) |> imap(x -> split!(fields, x, UInt8('\t'))))
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
    refseq::LongDNASeq
    asmseq::LongDNASeq
    insignificant::BitVector
    aln::PairwiseAlignment{LongDNASeq, LongDNASeq}
    orfdata::Vector{ORFData}
end

function ReferenceAssembly(ref::Reference, asm::Assembly)
    aln = pairalign(OverlapAlignment(), seq, ref, ALN_MODEL).aln
    orfdata = map(reference.proteins) do protein
        ORFData(protein, aln, reference)
    end
    ReferenceAssembly(ref.seq, asm.seq, asm.insignificant, aln, orfdata)
end

# This probably needs to be extended by letting "extras" be some new type with
# customizable behaviour

#=
function segment_report_lines(refasm::ReferenceAssembly,
    errormsgs::Dict{String, Vector{ErrorMessage}},
)::Tuple{Bool, Vector{String}}
    is_ok = true
    lines = String[]

    # Check for presence
    if is_error(maybe_assembly) || is_error(maybe_reference)
        push!(lines, "ERROR: Missing segment")
        return (false, lines)
    end
    assembly, reference = unwrap(maybe_assembly), unwrap(maybe_reference)

    # Check segment length and return early if way too short
    if length(assembly.seq) < 2 * TERMINAL + 1
        push!(lines, beginning * "\n\t\tERROR: Sequence or depths length: $(length(assembly.seq))")
        return false
    end
=#




