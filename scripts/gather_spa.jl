# In this script, it gathers all SPA hits, then creates new indexable
# FASTA files for each basename based on what segments is SPA aligns to

using ErrorTypes
using FASTX
using BioSequences

@enum Segment::UInt8 PB1 PB2 PA HA NP NA MP NS
const _STR_SEGMENT = Dict(map(i -> string(i)=>i, instances(Segment)))
Segment(s::AbstractString) = get(() -> error("Unknown segment: \"$(s)\""), _STR_SEGMENT, s)

# Gets the first column of the second line, or None if there is only the header
# (i.e. no good matches)
function readspa(path::AbstractString)::Option{UInt}
    open(path) do file
        lines = eachline(file)
        (header, _) = iterate(lines) # header
        @assert startswith(header, "#Template\tNum\t")
        nextline = iterate(lines)
        nextline === nothing && return none
        line = strip(first(nextline))
        isempty(line) && return none
        num = split(line, '\t')[2]
        Thing(parse(UInt, num))
    end
end

"Get a dict: basename => (HA => 5) for all present, mapping segments"
function make_numdict(alndir::AbstractString)
    numdict = Dict{String, Vector{Tuple{Segment, UInt}}}()
    for basename in readdir(alndir)
        v = Tuple{Segment, UInt}[]
        for segment in instances(Segment)
            num = @unwrap_or readspa(joinpath(alndir, basename, string(segment) * ".spa")) continue
            push!(v, (segment, num))
        end
        numdict[basename] = v
    end
    numdict
end

function collect_sequences(refdir::AbstractString, numdict::Dict{String, Vector{Tuple{Segment, UInt}}})
    # First get a Dict("HA" => Set(1, 2, 3)) of present nums
    present = Dict{Segment, Set{UInt}}()
    for (basename, vec) in numdict
        for (segment, num) in vec
            push!(get!(Set{UInt}, present, segment), num)
        end
    end

    # Then iterate over all segments, fetching the ones in the set
    records = Dict{Segment, Dict{UInt, FASTA.Record}}()
    for (segment, set) in present
        dict = Dict{UInt, FASTA.Record}()
        records[segment] = dict
        fastapath = joinpath(refdir, string(segment) * ".fna")
        open(FASTA.Reader, fastapath) do reader
            isempty(set) && return nothing
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
            header = FASTA.identifier(record) * '_' * string(segment)
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
    main(ARGS[1], ARGS[2])
end
