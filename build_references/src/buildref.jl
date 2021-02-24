module t

#=
TODO:

* Make sure CD-HIT is installed

=#

using Pkg
Pkg.activate(".")
#Pkg.instantiate()

using BioSequences
using CodecZlib
using ErrorTypes
using FASTX
using Serialization
using Transducers
using UnicodePlots

using Base: RefValue

"Complete the entire workflow. Get an overview of all that happens here."
function main()
    # Get the data from the internet
    donwload_influenza_data()

    # Clean the data
    isdir("results") || mkdir("results")
    clean_the_data(joinpath("results", "genomeset.clean.dat"), joinpath("raw", "genomeset.dat.gz"))

    # Parse genomeset to SegmentData with all basic info
    segment_data = open(parse_cleaned_genomeset, joinpath("results", "genomeset.clean.dat"))
    
    # Load influenza.fna and update the seq field of the SegmentData
    open(joinpath("raw", "influenza.fna.gz")) do io
        add_sequences!(segment_data, GzipDecompressorStream(io))
    end

    # Load influenza.dat and update the ORF fields of SegmentData

    # Series of filters on the flumeta

    # Add in the human ones for good measure

    # Serialize to files
end

function donwload_influenza_data(force=false)
    isdir("raw") && !force && return nothing

    mkdir("raw")
    ftp_address = "https://ftp.ncbi.nih.gov/genomes/INFLUENZA/"
    for filename in [
        "genomeset.dat.gz",
        "influenza.dat.gz",
        "influenza.fna.gz",
        ]
        download(joinpath(ftp_address, filename), joinpath("raw", filename))
    end
end

function clean_the_data(outpath, inpath)
    file = open(inpath)
    outfile = open(outpath, "w")
    decompressed = GzipDecompressorStream(file)

    for line in eachline(decompressed) |> Map(strip) ⨟ Filter(!isempty)
        fields = split(line, '\t')
        
        # Check correct number of fields
        if length(fields) != 11
            println("\"$line\"")
            error()
        end

        gi, host, segment, subtype, country, year, len, name, age, gender, group = fields

        # Filter subtype
        subtype_upper = uppercase(subtype)
        subtype = if startswith(subtype_upper, "MIXED")
            ""
        elseif subtype in ["H1", "H11N9/N2", "H1N2v", "H3N2v"]
            ""
        elseif subtype == "H3N6,H3"
            "H3N6"
        elseif subtype == "H6N1,H6"
            "H6N1"
        else
            subtype_upper
        end

        # Filter name
        nn = parse_name(name)
        name = if is_error(nn)
            ""
        else
            unwrap(nn)
        end

        # Filter year
        year = if year in ["NON", "Unknown", "unknown"]
            ""
        else
            year
        end 

        println(outfile, join([gi, host, segment, subtype, country, year, len, name, age, gender, group], '\t'))
    end

    close(decompressed)
    close(outfile)
end

"Parses a name in genomeset.dat during cleaning, returning none if it's malformed.
If it's unexpectedly malformed, throw an error"
function parse_name(s::Union{String, SubString})::Option{String}
    isempty(s) && return none

    occursin(r"^Influenza [AB] [vV]irus", s) || error(s)
    ncodeunits(s) == 17 && return none
    rest = SubString(s, 18 + (codeunit(s, 18) == UInt(' ')):ncodeunits(s))

    isempty(rest) && return none
    rest2 = if codeunit(rest, 1) == UInt8('(') && codeunit(rest, ncodeunits(rest)) == UInt8(')')
        SubString(rest, 2:ncodeunits(rest)-1)
    else
        rest
    end

    m = match(r"\([Hh]\d+[nN]\d+\)$", rest2)
    rest3 = if m === nothing
        rest2
    else
        SubString(rest2, 1:ncodeunits(rest2) - ncodeunits(m.match))
    end
    
    return some(String(rest3))
end

@enum Segment::UInt8 PB1 PB2 PA HA NP NA MP NS

function parse_from_integer(::Type{Segment}, s::AbstractString)::Option{Segment}
    y = tryparse(UInt8, s)
    y === nothing && return none
    (iszero(y) | (y > 0x08)) && return none
    some(reinterpret(Segment, y - 0x01))
end

const SEGMENT_STR_DICT = Dict(string(s)=>s for s in instances(Segment))
function parse_from_name(::Type{Segment}, s::AbstractString)::Option{Segment}
    y = get(SEGMENT_STR_DICT, s, nothing)
    y === nothing && return none
    some(y)
end

@enum Species::UInt8 human swine avian other

function Species(s::Union{String, SubString})::Species
    s == "Human" && return human
    s == "Swine" && return swine
    s == "Avian" && return avian
    return other
end

struct SubType
    H::UInt8
    N::UInt8
end

Base.show(io::IO, s::SubType) = print(io, 'H', s.H, 'N', s.N)

@enum FluSpecies::UInt8 InfluenzaA InfluenzaB

@enum ProteinVariant::UInt8 begin
    pb2
    pb1
    n40 # very rarely seen
    pb1f2 # inf A only
    pa
    pax # inf A only
    ha
    np
    na
    nb # inf B only
    m1
    m2 # inf A only
    bm2 # inf B only
    ns1
    nep
end

struct Protein
    variant::ProteinVariant
    orfs::Vector{UnitRange{UInt16}}
end

struct SegmentData
    gi::String
    host::Species
    fluspecies::FluSpecies
    segment::Segment
    subtype::SubType
    year::Int
    len::Int
    name::String
    proteins::Vector{Protein}
    seq::RefValue{Option{LongDNASeq}}
end

function parse_subtype(s::Union{String, SubString})::Option{SubType}
    isempty(s) && return none
    m = match(r"^H(\d+)N(\d+)$", s)
    m === nothing && error("UNKNOWN SUBTYPE: $s")
    some(SubType(parse(UInt8, m[1]), parse(UInt8, m[2])))
end

function parse_cleaned_genomeset(io::IO)::Vector{SegmentData}
    result = SegmentData[]
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        data = @unwrap_or parse(SegmentData, line) continue
        push!(result, data)
    end
    result
end

function Base.parse(::Type{SegmentData}, line::Union{String, SubString{String}})::Option{SegmentData}
    fields = split(line, '\t')
    subtype = @? parse_subtype(fields[4])
    year = @? parse_year(fields[6])
    host = Species(fields[2])
    fluspecies = @? parse_fluspecies(fields[8])
    gi = String(fields[1])
    name = String(fields[8])
    segment = @? parse_from_integer(Segment, fields[3])
    len = parse(Int, fields[7])
    some(SegmentData(gi, host, fluspecies, segment, subtype, year, len, name, Protein[], Ref(none(LongDNASeq))))
end

function parse_year(s::Union{String, SubString})::Option{Int}
    isempty(s) && return none
    pos_slash_found = findfirst(isequal('/'), s)
    last_byte = pos_slash_found === nothing ? ncodeunits(s) : pos_slash_found - 1
    return some(parse(Int, SubString(s, 1:last_byte)))
end

function parse_fluspecies(name::Union{String, SubString{String}})::Option{FluSpecies}
    if startswith(name, "A/")
        return some(InfluenzaA)
    elseif startswith(name, "B/")
        return some(InfluenzaB)
    else
        return none
    end
end

function add_sequences!(segment_data::Vector{SegmentData}, io::IO)
    gb_map = Dict(data.gi => data for data in segment_data)
    n_updated = 0
    @assert length(segment_data) == length(gb_map) "GB identifers not unique"
    for record in FASTA.Reader(io)
        # headers begins with e.g. gi|59292|gb|X53029|Influenza A virus
        gb_accession = split(FASTA.identifier(record)::String, '|')[4]
        haskey(gb_map, gb_accession) || continue
        data = gb_map[gb_accession]
        @assert is_error(data.seq[])
        data.seq[] = some(FASTA.sequence(LongDNASeq, record))
        n_updated += 1
    end
    println("Updated $n_updated/$(length(gb_map)) of records")
    segment_data
end

    

isinteractive() || main()

end # module