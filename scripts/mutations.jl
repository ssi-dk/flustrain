#=
TODO:
    Baloxavir

=#
module t

using Influenza
using SequenceVariation
using ErrorTypes
using BioSymbols
using BioSequences
using FASTX

const ProteinTuple{T} = NTuple{length(instances(Protein)), T}

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

function load_diffs(io::IO)
    eachline(io) |> imap(strip) |> ifilter(!isempty) |> imap() do line
        parse(Diff{LongAminoAcidSeq, AminoAcid}, line)
    end |> collect
end

function load_seq(path::AbstractString)
    record = only(open(collect, FASTA.Reader, path))
    FASTA.sequence(LongAminoAcidSeq, record)
end

function load_consensus(conspath::AbstractString
)::Dict{String, ProteinTuple{Option{LongAminoAcidSeq}}}
    result = Dict{String, ProteinTuple{Option{LongAminoAcidSeq}}}()
    record = FASTA.Record()
    vec = fill(none(LongAminoAcidSeq), length(instances(Protein)))
    for basename in readdir(conspath)
        fill!(vec, none(LongAminoAcidSeq))
        result[basename] = open(FASTA.Reader, joinpath(conspath, basename, "curated.faa")) do reader
            while !eof(reader)
                read!(reader, record)
                protein_str = split(FASTA.header(record)::String, '_')[end]
                protein = unwrap(parse(Protein, protein_str))
                protein_index = reinterpret(UInt8, protein) + 0x01
                seq::LongAminoAcidSeq = FASTA.sequence(LongAminoAcidSeq, record)
                vec[protein_index] = some(seq)
            end
            Tuple(vec)
        end
    end
    result
end


end # t
