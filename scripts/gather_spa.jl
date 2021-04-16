# In this script, it gathers all SPA hits, then creates new indexable
# FASTA files for each basename based on what segments is SPA aligns to

using InfluenzaCore
using ErrorTypes
using FASTX
using BioSequences

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

# Gets the first column of the second line, or None if there is only the header
# (i.e. no good matches)
function readspa(path::AbstractString)::Option{UInt}
    open(path) do file
        lines = eachline(file) |> imap(strip) |> ifilter(!isempty) |> imap(x -> split(x, '\t'))
        fields = zip(lines, 1:2) |> imap(first) |> collect
        @assert fields[1][1:2] == ["#Template", "Num"]
        length(fields) < 2 ? none : some(parse(UInt, fields[2][2]))
    end
end

"Get a dict: basename => (HA => 5) for all present, mapping segments"
function make_numdict(alndir::AbstractString)::Dict{String, Vector{Tuple{Segment, UInt}}}
    readdir(alndir) |> imap() do basename
        # Iterator of (segment, maybe_num)
        instances(Segment) |> imap() do segment
            (segment, readspa(joinpath(alndir, basename, string(segment) * ".spa")))
        end |>
        # Remove the nones (interestingly, Julia PR 38905 should let the compiler
        # remove the branch in the following unwrap due to this ifilter here. Cool stuff.)
        ifilter() do (segment, maybe_num)
            !is_error(maybe_num)
        end |>
        # Unwrap the nones
        imap() do (segment, maybe_num)
            (segment, unwrap(maybe_num))
        end |>
        # Convert iterator to (basename, [elements])
        (it -> (basename, collect(it)))
    end |>
    Dict
end

function collect_sequences(refdir::AbstractString, numdict::Dict{String, Vector{Tuple{Segment, UInt}}})
    # First get a Dict("HA" => Set(1, 2, 3)) of present nums
    present = Dict{Segment, Set{UInt}}()
    for (basename, vec) in numdict, (segment, num) in vec
        push!(get!(Set{UInt}, present, segment), num)
    end

    # Then iterate over all segments, fetching the ones in the set
    records = Dict{Segment, Dict{UInt, FASTA.Record}}()
    for (segment, set) in present
        isempty(set) && continue
        dict = Dict{UInt, FASTA.Record}()
        records[segment] = dict
        fastapath = joinpath(refdir, string(segment) * ".fna")
        open(FASTA.Reader, fastapath) do reader
            seqnum = UInt(0)
            record = FASTA.Record()
            lastnum = maximum(set)
            while !eof(reader)
                read!(reader, record)
                seqnum += UInt(1)
                if in(seqnum, set)
                    dict[seqnum] = copy(record)
                end
                seqnum == lastnum && break
            end
            if seqnum < lastnum
                error("FASTA $fastapath requested num $lastnum, has only $seqnum records")
            end
        end
    end
    return records
end

function dump_sequences(alndir::AbstractString, numdict::Dict{String, Vector{Tuple{Segment, UInt}}},
records::Dict{Segment, Dict{UInt, FASTA.Record}})
    for (basename, vec) in numdict
        writer = open(FASTA.Writer, joinpath(alndir, basename, "cat.fna"))
        for (segment, num) in vec
            record = records[segment][num]
            header = FASTA.identifier(record)::String * '_' * string(segment)
            newrecord = FASTA.Record(header, FASTA.sequence(LongDNASeq, record))
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
    if length(ARGS) != 2
        error("Usage: julia gather_spa.jl alndir refdir")
    end
    main(ARGS[1], ARGS[2])
end
