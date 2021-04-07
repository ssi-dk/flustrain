USAGE = """
update.jl - update influenza references

Purpose: Add new sequences to the collection of references used in the pipeline

\$ julia update.jl fluoutput species year subtypes passed
Where "fluoutput" is the top-level directory of a flu pipeline output
"species" is human, avian or swine,
subtypes is the path of the serialized "subtypes.jls" file from the build_references
passed is a file with "basename\\ssegment" lines of passed segments.

Writes a FASTA file with the proper name out
"""

using BioSequences
using FASTX
using InfluenzaCore
using Serialization
using ErrorTypes

# Format is ID|SPECIES|SEGMENT|SUBTYPE|YEAR|NAME|CLADE

function main(dir::AbstractString,
    speciesstr::AbstractString,
    yearstr::AbstractString,
    subtypepath::AbstractString,
    passedstr::AbstractString
)
    year = parse(Int, yearstr)
    1850 < year < 2050 || error("Year must be in 1850:2050")
    species = get_species(speciesstr)
    passed = parse_passed(passedstr)
    subdirs = readdir(joinpath(dir, "consensus"))
    if !issubset(keys(passed), subdirs)
        bad = first(setdiff(subdirs), keys(passed))
        error("Basename not found in consensus dir: $bad")
    end

    subtypedict = load_subtypes(subtypepath, dir, passed)

    # Read sequences
    records = FASTA.Record[]
    for (basename, segments) in passed
        open(joinpath(dir, "consensus", basename, "consensus.fna")) do io
            add_segments!(io, species, year, records, segments, subtypedict[basename])
        end
    end

    for record in records
        println(record)
    end
end

function get_species(s::AbstractString)
    s = lowercase(strip(s))
    if !in(s, ["human", "avian", "swine"])
        error("Species must be human, avian or swine")
    end
    s
end

function parse_passed(filename::AbstractString)
    open(filename) do file
        result = Dict{String, Vector{Segment}}()
        for line in eachline(file)
            stripped = strip(line)
            isempty(stripped) && continue
            
            basename, segmentstr, rest... = split(stripped)
            @assert isempty(rest)
            segment = expect(parse(Segment, segmentstr), "Error: Could not parse segment $segmentstr")
            push!(get!(Vector{Segment}, result, basename), segment)
        end
        result
    end
end

function load_subtypes(
    subtypepath::AbstractString,
    topdir::AbstractString,
    passed::Dict{String, Vector{Segment}},
)::Dict{String, Vector{SubType}}
    subtypedict = Dict{String, SubType}(open(deserialize, subtypepath))
    res = Dict{String, Vector{SubType}}()
    for (basename, segments) in passed
        headers = filter(collect(eachline(joinpath(topdir, "aln", basename, "kma1.fsa")))) do line
            startswith(line, '>')
        end
        pairs = map(headers) do header
            p = findlast(isequal('_'), header)
            refname = header[2:p-1]
            segmentstr = header[p+1:end]
            segment = unwrap(parse(Segment, segmentstr))
            (refname, segment)
        end
        subtypes = map(segments) do segment
            index = findfirst(i -> last(i) == segment, pairs)
            subtypedict[first(pairs[index])]
        end
        res[basename] = subtypes
    end
    res
end

function add_segments!(
    io::IO,
    species::String,
    year::Integer,
    records::Vector{FASTA.Record},
    segments::Vector{Segment},
    subtypes::Vector{SubType}
)
    reader = FASTA.Reader(io)
    n_seqs = 0
    for record in reader
        id = FASTA.identifier(record)
        basename, segmentstr = rsplit(id, '_', limit=2)
        segment = unwrap(parse(Segment, segmentstr))
        if segment in segments
            subtype = subtypes[findfirst(==(segment), segments)]
            subtype_str = if subtype.data isa SubTypes.InfluenzaA
                "H$(subtype.data._1)N$(subtype.data._2)"
            else
                string(typeof(subtype.data).name.name)
            end
            new_header = "$(basename)_$segment|$species|$segment|$subtype_str|$year|$basename|"
            new_record = FASTA.Record(new_header, FASTA.sequence(LongDNASeq, record))
            push!(records, new_record)
            n_seqs += 1
        end
    end
    @assert n_seqs == length(segments)
    records
end

# fluoutput species year subtypes passed
if length(ARGS) != 5
    println(USAGE)
    error()
end
main(ARGS...)