# This script is used to re-build the ref directory.
# To use this, install cd-hit, and get the influenza database from NCBI.
# They have an FTP link. I did this 2021-01-21

#=
Accession (string) => FASTA
Accession (string) => Protein

=#

using FASTX
using BioSequences
using ErrorTypes
using Transducers
using Serialization

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
    open("/Users/jakobnissen/Downloads/INFLUENZA/processed/genomeset.filt.dat", "w") do outfile
        clean_data(infile, outfile)
    end
end

###########################

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

const SEGMENT_STR_DICT = Dict(string(s)=>s for s in instances(Segment))

function parse_from_integer(::Type{Segment}, s::AbstractString)::Option{Segment}
    y = tryparse(UInt8, s)
    y === nothing && return none
    (iszero(y) | (y > 0x08)) && return none
    Thing(reinterpret(Segment, y - 0x01))
end

function parse_from_name(::Type{Segment}, s::AbstractString)::Option{Segment}
    y = get(SEGMENT_STR_DICT, s, nothing)
    y === nothing && return none
    Thing(y)
end

@enum ProteinVariant::UInt8 begin
    pb2
    pb1
    pb1fa
    pa
    pax
    ha
    np
    na
    m1
    m2
    ns1
    nep
end

# Approx +/- 100 bp for each ORF. For detecting which ORF is which
const ORF_LEN_RANGE = UnitRange{UInt16}[
    2100:2325,
    2100:2325,
    6:1000,
    2050:2300,
    600:800,
    1600:1800,
    1400:1600,
    1300:1500,
    650:850,
    200:400,
    550:750,
    250:450
]

orf_len_range(v::ProteinVariant) = @inbounds ORF_LEN_RANGE[reinterpret(UInt8, v) + 0x01]

const VARIANTS = [
    [pb2],
    [pb1, pb1fa],
    [pa, pax],
    [ha],
    [np],
    [na],
    [m1, m2],
    [ns1, nep]
]

variants(s::Segment) = @inbounds VARIANTS[reinterpret(UInt8, s) + 0x01]

# No known influenza proteins has more than two ORFs
struct Protein
    variant::ProteinVariant
    orfs::Union{Tuple{UnitRange{UInt16}}, NTuple{2, UnitRange{UInt16}}}
end

# No infleunza segment has more than two known proteins
struct Accession
    name::String
    segment::Segment
    proteins::Union{Tuple{Protein}, NTuple{2, Protein}}
end

#########################################################
"detects the ProteinVariant based on length of ORF"
function infer_proteinvariant(segment::Segment, orflen::Integer)::Option{ProteinVariant}
    for var in variants(segment)
        if orflen in orf_len_range(var)
            return Thing(var)
        end
    end
    return none
end

#########################################################

#### LOAD IN ORF DATA from "influenza.dat"
# Looks like one of these lines:
# (gb|AB000613:4-731, 960)
# gb|AB000721:710-1129
# gb|AB266385:<1-755
function parse_range(s::AbstractString)::Option{UnitRange{UInt16}}
    p2 = findfirst(isequal('-'), s)
    # If it's not a range, it must be a single number
    if p2 === nothing
        n = parse(UInt16, s)
        Thing(n:n)
    else
        start = parse(UInt16, @view s[1:p2-1])
        stop = parse(UInt16, s[p2+1:ncodeunits(s)])
        Thing(start:stop)
    end
end

# This returns the ORFS. When checking length, we can verify which of the ORFS it is.
function parse_orf(s::AbstractString, segment::Segment,
)::Option{Protein}
    # If it looks like gb|AB266090:<411->632, the ORF is not present
    # in the reference, and we skip it
    if occursin('>', s) || occursin('<', s)
        return none
    end

    # Multiple ORFs in a gene makes them enclosed in brackets
    s1 = strip(s, ['(', ')'])
    p = findfirst(isequal(':'), s1)
    p === nothing && return none
    s2 = @view s1[p+1:ncodeunits(s1)]
    (range, len) = if occursin(',', s2)
        orfstrings = split(s2, ", ")
        length(orfstrings) == 2 || return none
        r1 = @? parse_range(orfstrings[1]) 
        r2 = @? parse_range(orfstrings[2])
        (r1, r2), (length(r1) + length(r2))
    else
        r1 = @? parse_range(s2)
        (r1,), length(r1)
    end
    # Return none if it's not divisible by three and thus cant be ORF
    iszero(len % 3) || return none
    var = @? infer_proteinvariant(segment, len)
    Thing(Protein(var, range))
end

# First col is the segment accession, fields 3, 5, 7 etc. are protein ORFS
function load_orf_line(v::Vector{Protein}, line::AbstractString, segment::Segment)::Option{String}
    fields = split(line, '\t')
    length(fields) < 3 && return none
    empty!(v)
    for orf_field in fields[3:2:end]
        push!(v, @?(parse_orf(orf_field, segment)))
    end
    return Thing(String(first(fields)))
end

function load_orf_data(io::IO, segmentof::Dict{String, Segment})
    orfs = Dict{String, Union{Tuple{Protein}, NTuple{2, Protein}}}()
    v = Vector{Protein}()
    for line in eachline(io)
        id = first(split(line, '\t'))
        segment = get(segmentof, id, nothing)
        segment === nothing && continue
        maybe_key = load_orf_line(v, line, segment)
        is_none(maybe_key) && continue
        haskey(orfs, unwrap(maybe_key)) && error("Duplicate key $(unwrap(maybe_key))!")
        orfs[unwrap(maybe_key)] = length(v) == 2 ? Tuple(v) : (first(v),)
    end
    orfs
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
    subtype::SubType
    year::Int
    len::Int
    name::String
end

function parse_flumeta(s::Union{String, SubString})::Option{FluMeta}
    fields = split(s, '\t')
    subtype = @? parse_subtype(fields[4])
    year = @? parse_year(fields[6])
    host = Species(fields[2])
    gi = String(fields[1])
    name = String(fields[8])
    segment = @? parse_from_integer(Segment, fields[3])
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

records = open(parse_all, "/Users/jakobnissen/Downloads/INFLUENZA/processed/genomeset.filt.dat")
segmentof = Dict(m.gi => m.segment for m in records)
orf_dict = open("/Users/jakobnissen/Downloads/INFLUENZA/raw/influenza.dat") do io
    load_orf_data(io, segmentof)
end

##################### A series of filters to filter the flumeta
###############################################################
const N_PROTEINS = [(1,), (1, 2), (2,), (1,), (1,), (1,), (2,), (2,)]
n_proteins(s::Segment) = @inbounds N_PROTEINS[reinterpret(UInt8, s) + 0x01]

filter!(records) do record
    proteins = get(orf_dict, record.gi, nothing)

    # Has to have parse-able ORFs
    proteins === nothing && return false
    n = isa(proteins, Tuple) ? length(proteins) : 1

    # Has to have the correct number of ORFs
    n in n_proteins(record.segment)
end

# Load in seqs by genbank acc number
by_id = open(FASTA.Reader, "/Users/jakobnissen/Downloads/INFLUENZA/raw/influenza.fna") do reader
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

########## Check lengths
bysegment = Dict(s => FluSeq[] for s in instances(Segment))
foreach(fluseqs) do fluseq
    push!(bysegment[fluseq.meta.segment], fluseq)
end

# Plot lengths to see outliers
for (k, v) in bysegment
    show(histogram(map(x -> length(x.seq), v), title=string(k)))
end

# This is manually made by looking at the plots above, in order
# to remote truncated seqs or seqs with junk in them
trimrange = Dict(
    NP => 1490:1570,
    HA => 1680:1780,
    MP => 980:1030,
    PB1 => 2270:2345,
    PA => 2150:2240,
    NA => 1350:1470,  # very wide distribution??
    NS => 820:895,
    PB2 => 2280:2345
)

filter!(fluseqs) do fluseq
    in(length(fluseq.seq), trimrange[fluseq.meta.segment])
end

#### NOW WE CAN DEDUPLICATE ET CETERA

# Create a unfilt for each segment
for species in [avian, swine]
    hostseqs = filter(x -> x.meta.host == species, fluseqs)
    for segment in instances(Segment)
        seqs = filter(hostseqs) do seq
            seq.meta.segment == segment
        end
        open("/Users/jakobnissen/Downloads/INFLUENZA/processed/$species/$segment.fna", "w") do file
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
        inpath = "/Users/jakobnissen/Downloads/INFLUENZA/processed/$species/$segment.fna"
        outpath = inpath * ".cdhit"
        command = `$exec -i $inpath -o $outpath -aS 0.9 -c 0.95`
        run(command)
    end
end

#### For each of the kept sequences, we get the ORF and serialize it to a file!.
for species in [avian, swine]
    for segment in instances(Segment)
    # Load the kept seqs!
        deduplicated_headers = open("/Users/jakobnissen/Downloads/INFLUENZA/processed/$species/$segment.fna.cdhit") do file
            [line[2:end] for line in eachline(file) if startswith(line, '>')]
        end

        T = Tuple{UInt8, Union{Tuple{UnitRange{UInt16}}, NTuple{2, UnitRange{UInt16}}}}
        kept_orfs = Vector{Tuple{String, Union{Tuple{T}, NTuple{2, T}}}}()
        for acc in deduplicated_headers
            proteins = orf_dict[acc]
            in(acc, kept_accessions) || continue
            v = map(proteins) do protein
                (reinterpret(UInt8, protein.variant), protein.orfs)
            end
            push!(kept_orfs, (acc, v))
        end

        open("/Users/jakobnissen/Downloads/INFLUENZA/processed/$species/$segment.orf.jls", "w") do file
            serialize(file, kept_orfs)
        end
    end
end

# Move to the correct directory (commented out for safety)

cd("/Users/jakobnissen/Downloads/INFLUENZA/")
target_dir = "/Users/jakobnissen/Documents/ssi/projects/flupipe/ref/"
for species in [avian, swine]
    subdir = joinpath(target_dir, string(species))
    isdir(subdir) || mkdir(subdir)
    for segment in instances(Segment)
        # Move copied file
        cp("processed/$species/$segment.fna.cdhit",
        "$subdir/$segment.fna",
        force=true)

        # Move serialized ORFs
        cp("processed/$species/$segment.orf.jls",
        "$subdir/$segment.orfs.jls",
        force=true)
    end
end
