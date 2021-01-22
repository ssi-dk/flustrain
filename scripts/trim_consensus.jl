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

function is_possibility(needle::NucleotideSeq, haystack::NucleotideSeq)
    @inbounds for i in eachindex(needle)
        iscompatible(needle[i], haystack[i]) || return false
    end
    return true
end

function slide!(primer::NucleotideSeq, seq::NucleotideSeq, minlen::Int)
    length(seq) < length(primer) && return 0
    primer = copy(primer)
    seq = copy(seq)[1:length(primer)]
    length(primer) == length(seq) || throw(ArgumentError("Must match in lengths"))
    for overlap in length(primer):-1:minlen
        is_possibility(primer, seq) && return overlap
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

function load_consensus(path::String)
    open(FASTA.Reader, path) do reader
        record, _ = iterate(reader)
        seq = FASTA.sequence(LongDNASeq, record)
        header = get_header(record)
        iterate(reader) === nothing || error("Multiple records in file $path")
        (header, seq)
    end
end

function remove_primers(seq::NucleotideSeq, primers::Vector{<:Tuple{String, NucleotideSeq}}, minlength::Int)
    overlaps = [slide!(p, seq, minlength) for (h,p) in primers]
    if maximum(overlaps) > 0
        i = argmax(overlaps)
        overlap = overlaps[i]
        header = primers[i][1]
        println(stdout, "Primer $header found with overlap $overlap")
        seq = seq[1+overlap:end]
    end
    return seq
end

function trim_consensus(primerpath::String, consensuspath::String, output::String, minlength::Int)
    primers = load_primers(primerpath)
    header, seq = load_consensus(consensuspath)

    if !isempty(primers)
        seq = remove_primers(seq, primers, minlength)
        seq = remove_primers(reverse_complement(seq), primers, minlength)
        seq = reverse_complement(seq)
    end

    open(FASTA.Writer, output) do writer
        write(writer, FASTA.Record(header, seq))
    end 
    
end 

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 4
        error("Usage: julia trim_consensus.jl primers.fna consensus.fna output.fna minlength")
    end
    minlength = parse(Int, ARGS[4])
    minlength < 1 && error("Minlength must be one")

    trim_consensus(ARGS[1], ARGS[2], ARGS[3], minlength)
end
