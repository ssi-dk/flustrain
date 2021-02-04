# This script scans a dir of KMA-assembled consensus seqs, and em out
# by their segment

using FASTX
using BioSequences

function main(alnpath::AbstractString, consensuspath::AbstractString, segments::AbstractString)
    segments = split(segments, ',')
    basenames = readdir(alnpath)
    for basename in basenames
        fastapath = joinpath(alnpath, basename, "kma2.fsa")
        records = open(collect, FASTA.Reader, fastapath)
        dict = Dict(FASTA.identifier(record) => FASTA.sequence(LongDNASeq, record) for record in records)
        if !isempty(setdiff(keys(dict), segments))
            badheader = first(setdiff(keys(dict), segments))
            error("Unknown segment header in FASTA $(fastapath): $badheader")
        end
        for segment in segments
            open(FASTA.Writer, joinpath(consensuspath, basename, segment * ".fna")) do writer
                if haskey(dict, segment)
                    newrecord = FASTA.Record(basename * '_' * segment, dict[segment])
                    write(writer, newrecord)
                end
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], ARGS[2], ARGS[3])
end
