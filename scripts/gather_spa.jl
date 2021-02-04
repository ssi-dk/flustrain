# In this script, it gathers all SPA hits, then creates new indexable
# FASTA files for each basename based on what segments is SPA aligns to

using ErrorTypes
using FASTX
using BioSequences

function readspa(path::AbstractString)::Option{UInt}
    open(path) do file
        lines = eachline(file)
        iterate(lines) # header
        nextline = iterate(lines)
        nextline === nothing && return none
        line = strip(first(nextline))
        isempty(line) && return none
        num = split(line, '\t')[2]
        Thing(parse(UInt, num))
    end
end

"Get a basename => [\"HA\" => 42, \"NA\" => 11] for all present, mapping segments"
function make_numdict(alndir::AbstractString)
    numdict = Dict{String, Vector{Tuple{String, UInt}}}()
    basenames = readdir(alndir)
    for basename in basenames
        spas = filter(readdir(joinpath(alndir, basename))) do filename
            endswith(filename, ".spa")
        end
        sort!(spas)
        segments = map(spas) do filename
            filename[1 : findfirst(isequal('.'), filename) - 1]
        end
        vec = Vector{Tuple{String, UInt}}()
        numdict[basename] = vec
        for (segment, spa) in zip(segments, spas)
            num = readspa(joinpath(alndir, basename, spa))
            is_none(num) && continue
            push!(vec, (segment, unwrap(num)))
        end
    end
    numdict
end

function collect_sequences(refdir::AbstractString, numdict::Dict{String, Vector{Tuple{String, UInt}}})
    # First get a Dict("HA" => Set(1, 2, 3)) of present nums
    present = Dict{String, Set{UInt}}()
    for (basename, vec) in numdict
        for (segment, num) in vec
            haskey(present, segment) || (present[segment] = Set{UInt}())
            push!(present[segment], num)
        end
    end

    # Then iterate over all segments, fetching the ones in the set
    records = Dict{String, Dict{UInt, FASTA.Record}}()
    for (segment, set) in present
        dict = Dict{UInt, FASTA.Record}()
        records[segment] = dict
        refpath = joinpath(refdir, segment * ".fna")
        open(FASTA.Reader, refpath) do reader
            seqnum = 0
            record = FASTA.Record()
            while !eof(reader)
                read!(reader, record)
                seqnum += 1
                if in(seqnum, set)
                    dict[seqnum] = copy(record)
                end
            end
        end
    end

    return records
end

function dump_sequences(alndir::AbstractString, numdict::Dict{String, Vector{Tuple{String, UInt}}},
records::Dict{String, Dict{UInt, FASTA.Record}})
    for (basename, vec) in numdict
        writer = open(FASTA.Writer, joinpath(alndir, basename, "cat.fna"))
        for (segment, num) in vec
            record = records[segment][num]
            newrecord = FASTA.Record(segment, FASTA.sequence(LongDNASeq, record))
            write(writer, newrecord)
        end
        close(writer)
    end
end

function main(alndir::AbstractString, refdir::AbstractString)
    numdict = make_numdict(alndir)
    records = collect_sequences(refdir, numdict)
    dump_sequences(alndir, numdict, records)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], ARGS[2])
end
