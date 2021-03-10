struct Assembly
    insignificant::BitVector
    seq::LongDNASeq
    accession::String
end

function load_assembly(assemblypath::AbstractString, kma::Bool)::SegmentTuple{Option{Assembly}}
    open(FASTA.Reader, assemblypath) do reader
        result = fill(none(Assembly), length(instances(Segment)))
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
            segment = unwrap(parse(Segment, v[2]))
            segment_index = reinterpret(UInt8, segment) + 0x01
            seq = FASTA.sequence(LongDNASeq, record)
            @assert length(record.sequence) == length(seq)
            if kma
                insignificant = BitVector([in(i, UInt8('a'):UInt8('z')) for i in @view record.data[record.sequence]])
            else
                insignificant = falses(length(record.sequence))
            end
            is_error(result[segment_index]) || error("Segment $segment present twice in $assemblypath")
            result[segment_index] = some(Assembly(insignificant, seq, accession))
        end
        SegmentTuple(result)
    end
end