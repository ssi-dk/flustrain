struct Assembly
    segment::Segment
    insignificant::BitVector
    seq::LongDNASeq
    accession::String
end

function load_assembly(assemblypath::AbstractString, kma::Bool)::SegmentTuple{Option{Assembly}}
    open(FASTA.Reader, assemblypath) do reader
        result = fill(none(Assembly), N_SEGMENTS)
        foreach(reader) do record
            id = FASTA.identifier(record)::String
            v = rsplit(id, '_', limit=2 + !kma)
            if kma
                length(v) == 2 || error("Found header \"$(id)\", expected pattern \"HEADER_SEGMENT\"")
            else
                if !(length(v) == 3 && last(v) == "segment0")
                    error("Found header \"$(id)\", expected pattern \"HEADER_SEGMENT_segment0\"")
                end
            end
            accession = first(v)
            segment::Segment = tryparse(Segment, v[2])
            segment_index = reinterpret(UInt8, segment) + 0x01
            seq = FASTA.sequence(LongDNASeq, record)
            @assert length(record.sequence) == length(seq)
            if kma
                insignificant = BitVector([in(i, UInt8('a'):UInt8('z')) for i in @view record.data[record.sequence]])
            else
                insignificant = falses(length(record.sequence))
            end
            is_error(result[segment_index]) || error("Segment $segment present twice in $assemblypath")
            result[segment_index] = some(Assembly(segment, insignificant, seq, accession))
        end
        SegmentTuple(result)
    end
end

# Here, we load template identities from kma2.res. If a template has <99.5% identity,
# it means the first assembly did not converge.
function load_kma_file(resfilename::AbstractString)::SegmentTuple{Option{Float64}}
    open(resfilename) do io
        fields = Vector{SubString{String}}(undef, 11)
        lines = eachline(io)
        header, _ = iterate(lines)::NTuple{2, Any}
        result = fill(none(Float64), N_SEGMENTS)
        @assert header == "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value"
        for fields in lines |> Map(strip) ⨟ Filter(!isempty) ⨟ Map(x -> split!(fields, x, UInt8('\t')))
            segment::Segment = tryparse(Segment, strip(last(rsplit(first(fields), '_', limit=2))))
            segment_index = reinterpret(UInt8, segment) + 0x01
            is_error(result[segment_index]) || error("Segment $(string(segment)) present twice in file $resfilename")
            result[segment_index] = some(parse(Float64, fields[5]) / 100)
        end
        SegmentTuple(result)
    end
end

function check_second_kma_file(resfilename::AbstractString)::SegmentTuple{Option{ErrorMessage}}
    result = fill(none(ErrorMessage), N_SEGMENTS)
    for (i, maybe_identity) in enumerate(load_kma_file(resfilename))
        identity = @unwrap_or maybe_identity continue
        if identity < 0.995
            percent = identity * 100
            result[i] = some(ErrorMessage(important, "Assembly not converged, identity at $(percent) %"))
        end
    end
    SegmentTuple(result)
end