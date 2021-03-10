struct ProteinORF
    var::Protein
    # a bitvector over its reference seq with 1s for coding nucleotides, incl last stop
    mask::BitVector
end

struct Reference
    segment::Segment
    seq::LongDNASeq
    proteins::Vector{ProteinORF}
end

function load_references(references::Set{Tuple{Segment, String}}, refdir::AbstractString
    )::Dict{Tuple{Segment, String}, Reference}
    # Group segments
    bysegment = Dict{Segment, Set{String}}()
    for (segment, accession) in references
        push!(get!(Set{String}, bysegment, segment), accession)
    end

    record = FASTA.Record()
    result = Dict{Tuple{Segment, String}, Reference}()
    for (segment, accessions) in bysegment
        seqpath = joinpath(refdir, "$segment.fna")
        orfpath = joinpath(refdir, "$segment.jls")

        # Get the sequence
        seqof = Dict{String, LongDNASeq}()
        open(FASTA.Reader, seqpath) do reader
            while !eof(reader)
                read!(reader, record)
                accession = FASTA.identifier(record)::String
                if in(accession, accessions)
                    seqof[accession] = FASTA.sequence(LongDNASeq, record)
                end
            end
        end

        # Get the proteins
        added_accessions = Set{String}()
        open(orfpath) do io
            for (accession, v) in deserialize(io)
                if in(accession, accessions)
                    if !haskey(seqof, accession)
                        error("Accession $accession missing from $seqpath")
                    end
                    seq = seqof[accession]
                    reference = Reference(segment, seq, ProteinORF[])
                    result[(segment, accession)] = reference
                    for (var_uint8, orf_tuple) in v
                        mask = falses(length(seq))
                        for orf in orf_tuple
                            @view(mask[orf]) .= true
                        end
                        push!(reference.proteins, ProteinORF(Protein(var_uint8), mask))
                    end
                    push!(added_accessions, accession)
                end
            end
        end
        
        # Check we added all accessions
        missing_acc = setdiff(accessions, added_accessions)
        if !isempty(missing_acc)
            error("Accession $(first(missing_acc)) missing from $orfpath")
        end
    end
    result
end

function load_references(
    maybe_asm_tuples::Vector{SegmentTuple{Option{Assembly}}},
    refdir::AbstractString
)::Vector{SegmentTuple{Option{Reference}}}
    # Create (asm => (segment, string)) dict
    accession_map::Dict{Assembly, Tuple{Segment, String}} = maybe_asm_tuples |>
    Cat() ⨟
    Filter(!is_error) ⨟
    Map() do maybe_asm
        asm = unwrap(maybe_asm)
        asm => (asm.segment, asm.accession)
    end |> Dict

    reference_map::Dict{Tuple{Segment, String}, Reference} = 
        load_references(Set(values(accession_map)), refdir)

    maybe_asm_tuples |> Map() do asm_tuple
        map(asm_tuple) do maybe_asm
            and_then(Reference, maybe_asm) do asm
                reference_map[(asm.segment, asm.accession)]
            end
        end
    end |> collect
end