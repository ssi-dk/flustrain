# Purpose: Trim primers off consensus seq

using FASTX
using BioSequences

function get_header(record::FASTA.Record)
    emptyid = isempty(record.identifier)
    emptydesc = isempty(record.description)
    if (emptyid & emptydesc)
        return ""
    elseif emptydesc
        return String(record.data[record.identifier])
    elseif emptyid
        return String(record.data[record.description])
    else
        return String(record.data[first(record.identifier):last(record.description)])
    end
end

function is_possibility(needle::NucleotideSeq, haystack::NucleotideSeq, maxerrs::Int)
    errs = 0
    @inbounds for i in eachindex(needle)
        if !iscompatible(needle[i], haystack[i])
            errs += 1
            errs > maxerrs && return false
        end
    end
    return true
end

function slide!(primer::NucleotideSeq, seq::NucleotideSeq, minlen::Int, fuzzylen::Int)
    length(seq) < length(primer) && return 0
    primer = copy(primer)
    seq = copy(seq)[1:length(primer)]
    length(primer) == length(seq) || throw(ArgumentError("Must match in lengths"))
    for overlap in length(primer):-1:minlen
        maxerrs = ifelse(overlap < fuzzylen, 0, 1)
        is_possibility(primer, seq, maxerrs) && return overlap
        popfirst!(primer)
        pop!(seq)
    end
    return 0
end

function load_primers(path::String)
    open(FASTA.Reader, path) do reader
        map(record -> (get_header(record), FASTA.sequence(LongDNASeq, record)), reader)
    end
end

function load_consensus(path::String)::Vector{Tuple{String, LongDNASeq}}
    recs = open(collect, FASTA.Reader, path)
    [(get_header(record), FASTA.sequence(LongDNASeq, record)) for record in recs]
end

function remove_primers(seqheader::String, seq::NucleotideSeq,
primers::Vector{<:Tuple{String, NucleotideSeq}}, minlength::Int, fuzzylen::Int
)
    overlaps = [slide!(p, seq, minlength, fuzzylen) for (h,p) in primers]
    if maximum(overlaps) > 0
        i = argmax(overlaps)
        overlap = overlaps[i]
        header = primers[i][1]
        println("Primer $header found in $seqheader with overlap $overlap")
        seq = seq[1+overlap:end]
    end
    return seq
end

function trim_consensus(primerpath::String, consensuspath::String, output::String,
minlength::Int, fuzzylen::Int)
    consensus = load_consensus(consensuspath)
    primers = load_primers(primerpath)
    open(FASTA.Writer, output) do writer
        for (header, seq) in consensus
            if !isempty(primers)
                seq = remove_primers(header, seq, primers, minlength, fuzzylen)
                seq = remove_primers(header, reverse_complement(seq), primers, minlength, fuzzylen)
                seq = reverse_complement(seq)
            end
            newrecord = FASTA.Record(header, seq)
            write(writer, newrecord)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 5
        error("Usage: julia trim_consensus.jl primers.fna consensus.fna output.fna minlength fuzzylen")
    end
    minlength = parse(Int, ARGS[4])
    minlength < 1 && error("Minlength must be one")
    fuzzylen = parse(Int, ARGS[5])
    fuzzylen >= minlength || error("Fuzzylen cannot be smaller than minlength")

    trim_consensus(ARGS[1], ARGS[2], ARGS[3], minlength, fuzzylen)
end
