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

function load_references(
    segment::Segment,
    refdir::AbstractString,
    accessions::Set{String},
)::Dict{String, Reference}
    seqpath = joinpath(refdir, "$segment.fna")
    orfpath = joinpath(refdir, "$segment.jls")
    result = Dict{String, Reference}()
    
    # Get the sequence
    seqof = Dict{String, LongDNASeq}()
    record = FASTA.Record()
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
                result[accession] = reference
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
    return result
end

function load_references(
    maybe_asm_tuples::Vector{SegmentTuple{Option{Assembly}}},
    refdir::AbstractString
)::Vector{SegmentTuple{Option{Reference}}}

    # Make a segment => [accessions ... ] dict
    accession_dict = Dict(s => Set{String}() for s in instances(Segment))
    foreach(maybe_asm_tuples) do maybe_asm_tuple
        for maybe_asm in maybe_asm_tuple
            asm = @unwrap_or maybe_asm continue
            push!(accession_dict[asm.segment], asm.accession)
        end
    end

    # Make a segment => Dict(accession => reference) nested dict
    reference_map = Dict(
        segment => load_references(segment, refdir, accessions)
        for (segment, accessions) in accession_dict
    )

    # Look up in the dict to return the result
    maybe_asm_tuples |> Map() do asm_tuple
        map(asm_tuple) do maybe_asm
            and_then(Reference, maybe_asm) do asm
                reference_map[asm.segment][asm.accession]
            end
        end
    end |> collect
end