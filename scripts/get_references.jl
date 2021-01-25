# This script is used to re-build the ref directory.
# To use this, install cd-hit, and get the influenza database from NCBI.
# They have an FTP link. I did this 2021-01-21

using FASTX
using BioSequences
using ErrorTypes
using Transducers


##### CLEAN THE DATA

format_error(s) =  error("Unknown format: \"$s\"")
function parse_name(s::Union{String, SubString})::Option{String}
    isempty(s) && return none

    occursin(r"^Influenza [AB] [vV]irus", s) || format_error(s)
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
    
    return Thing(String(rest3))
end

function clean_data(io::IO, out::IO)
    for line in (eachline(io) |> Filter(!isempty))
        fields = split(line, '\t')
        if length(fields) != 11
            println("\"$line\"")
            error()
        end
        gi, host, segment, subtype, country, year, len, name, age, gender, group = fields

        # Filter segment
        if !in(parse(Int, segment), 1:8)
            segment = ""
        end

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
        name = if is_none(nn)
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

        println(out, join([gi, host, segment, subtype, country, year, len, name, age, gender, group], '\t'))
    end
end

open("/Users/jakobnissen/Downloads/INFLUENZA/raw/genomeset.dat") do infile
    open("/Users/jakobnissen/Downloads/INFLUENZA/genomeset.filt.dat", "w") do outfile
        clean_data(infile, outfile)
    end
end

#### PARSE THE DATA

@enum Species::UInt8 begin
    human
    swine
    avian
    other
end

function Species(s::Union{String, SubString})::Species
    s == "Human" && return human
    s == "Swine" && return swine
    s == "Avian" && return avian
    return other
end

@enum Segment::UInt8 begin
    PB2
    PB1
    PA
    HA
    NP
    NA
    MP
    NS
end

Segment(s::Union{String, SubString}) = Segment(parse(UInt8, s) - 0x01)

struct SubType
    H::UInt8
    N::UInt8
end

Base.show(io::IO, s::SubType) = print(io, 'H', s.H, 'N', s.N)

function parse_subtype(s::Union{String, SubString})::Option{SubType}
    isempty(s) && return none
    m = match(r"^H(\d+)N(\d+)$", s)
    m === nothing && error("UNKNOWN SUBTYPE: $s")
    Thing(SubType(parse(UInt8, m[1]), parse(UInt8, m[2])))
end

function parse_year(s::Union{String, SubString})::Option{Int}
    isempty(s) && return none
    pos_slash_found = findfirst(isequal('/'), s)
    last_byte = pos_slash_found === nothing ? ncodeunits(s) : pos_slash_found - 1
    return Thing(parse(Int, SubString(s, 1:last_byte)))
end

struct FluMeta
    gi::String
    host::Species
    segment::Segment
    subtype::Option{SubType}
    year::Int
    len::Int
    name::String
end

function parse_flumeta(s::Union{String, SubString})::Option{FluMeta}
    fields = split(s, '\t')
    subtype = parse_subtype(fields[4])
    year = @? parse_year(fields[6])
    host = Species(fields[2])
    gi = String(fields[1])
    name = String(fields[8])
    segment = Segment(fields[3])
    len = parse(Int, fields[7])

    return Thing(FluMeta(gi, host, segment, subtype, year, len, name))
end

parse_all(io::IO) = parse_all(eachline(io) |> Map(strip) |> Filter(!isempty))

function parse_all(lines)
    metas = FluMeta[]
    for (i, line) in enumerate(lines)
        opt = parse_flumeta(line)
        is_none(opt) || push!(metas, unwrap(opt))
    end
    metas
end

records = parse_all(eachline("processed/genomeset.filt.dat"))

#### NOW WE CAN DEDUPLICATE ET CETERA
# The lengths are fine enough

# Load in seqs by genbank acc number
by_id = open(FASTA.Reader, "raw/influenza.fna") do reader
    record = FASTA.Record()
    records = Dict{String, LongDNASeq}()
    while !eof(reader)
        read!(reader, record)
        gi_id = String(split(FASTA.identifier(record), '|')[4])
        records[gi_id] = FASTA.sequence(LongDNASeq, record)
    end
    records
end

# Match the seq with metadata
struct FluSeq
    seq::LongDNASeq
    meta::FluMeta
end

fluseqs = FluSeq[]
for record in records # all can be matched
    seq = get(by_id, record.gi, nothing)
    seq === nothing && continue

    # Filter for Ns
    count(isambiguous, seq) > 4 && continue

    push!(fluseqs, FluSeq(seq, record))
end

# Create a unfilt for each segment
for species in [avian, swine]
    hostseqs = filter(x -> x.meta.host == species, fluseqs)
    for segment in instances(Segment)
        seqs = filter(hostseqs) do seq
            seq.meta.segment == segment
        end
        open("processed/$species/$segment.fna", "w") do file
            for seq in seqs
                println(file, '>', seq.meta.gi, '\n', seq.seq)
            end
        end
    end
end

# Now cluster for each. 95% identity.
exec = "/Users/jakobnissen/miniconda3/envs/bioinfo/bin/cd-hit-est"
for species in [avian, swine]
    Threads.@threads for segment in instances(Segment)
        inpath = "processed/$species/$segment.fna"
        outpath = inpath * ".cdhit"
        command = `$exec -i $inpath -o $outpath -aS 0.9 -c 0.95`
        run(command)
    end
end

for species in [avian, swine]
    for segment in instances(Segment)
        cp("processed/$species/$segment.fna.cdhit",
        "/Users/jakobnissen/Documents/ssi/projects/flupipe/ref/$species/$segment.fna",
        force=true)
    end
end