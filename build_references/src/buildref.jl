
using Pkg
Pkg.activate(".")
#Pkg.instantiate()

using BioSequences # commit 252720879620716f5f33a54e53ffcd9716706636
using CodecZlib
using ErrorTypes
using FASTX # commit e9ff73b65e6638ed740b2fce5a01b254ecb24c55
using Serialization
using Transducers
using UnicodePlots

using Base: RefValue

"Complete the entire workflow. Get an overview of all that happens here."
function main()
    # Get the data from the internet (except if it already exists)
    download_influenza_data(force=false)

    # Clean the data
    isdir("results") || mkdir("results")
    clean_the_data(joinpath("results", "genomeset.clean.dat"), joinpath("download", "genomeset.dat.gz"))

    # Parse genomeset to SegmentData with all basic info
    segment_data = open(parse_cleaned_genomeset, joinpath("results", "genomeset.clean.dat"))
    
    # Load influenza.fna and update the seq field of the SegmentData
    open(joinpath("download", "influenza.fna.gz")) do io
        add_sequences!(segment_data, GzipDecompressorStream(io))
    end

    # Load influenza_aa.dat.gz to get information on what kind of protein
    # the different ORFs make
    accession_protein_map = open(joinpath("download", "influenza_aa.dat.gz")) do io
        parse_inf_aa(GzipDecompressorStream(io))
    end

    # Load influenza.dat and update the ORF fields of SegmentData
    open(joinpath("download", "influenza.dat.gz")) do io
        add_orfs!(segment_data, accession_protein_map, GzipDecompressorStream(io))
    end

    # We remove the human sequences, because we want to add in a set of manually
    # curated human sequences. Also remove all species in we don't need.
    filter!(pair -> last(pair).host ∈ (avian, swine), segment_data)

    # Add in the human ones - these are annotated by The NCBI Influenza Virus Sequence Annotation Tool
    open(joinpath("raw", "annotation.tbl.gz")) do io
        add_human_sequences!(segment_data, GzipDecompressorStream(io), joinpath("raw", "human.fna.gz"))
    end

    # Series of filters on the SegmentData
    filter_segment_data!(segment_data)

    # Cluster using CD-hit for 95% identity
    # This is a Dict{Species, Dict{Segment, Set{String}}}, where String == accessions
    # (for human sequences, it's the "name" field instead)
    deduplicated = cd_hit_deduplicate(segment_data)

    # Serialize to files
    serialize_segments(segment_data, deduplicated)
end

function download_influenza_data(;force=false)
    ftp_address = "https://ftp.ncbi.nih.gov/genomes/INFLUENZA/"
    for filename in [
        "genomeset.dat.gz",
        "influenza.dat.gz",
        "influenza.fna.gz",
        "influenza_aa.dat.gz"
        ]
        download(joinpath(ftp_address, filename), joinpath("download", filename))
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
    B::UInt8 # 0 for inf A, 1 for victoria, 2 for yamagata
end

function Base.show(io::IO, s::SubType)
    iszero(s.B) && return print(io, 'H', s.H, 'N', s.N)
    print(io, s.B === 0x01 ? "Victoria" : "Yamagata")
end

@enum FluSpecies::UInt8 InfluenzaA InfluenzaB

# Note: Some of these variants may be transient names given in the current literature.
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
    m42 # uncommon
    ns1
    nep
    ns3 # uncommon
end

const _STR_PROTEINVARIANT = Dict(
    "PB2" => pb2,
    "PB1" => pb1,
    "N40" => n40,
    "PB1-F2" => pb1f2,
    "PA" => pa,
    "PA-X" => pax,
    "HA" => ha,
    "NP" => np,
    "NA" => na,
    "NB" => nb,
    "M1" => m1,
    "M2" => m2,
    "BM2" => bm2,
    "M42" => m42,
    "NS1" => ns1,
    "NS2" => nep,
    "NEP" => nep, # yeah, it has two names
    "NS3" => ns3,
)

# Parse how it's presented in the NCBI files
function Base.parse(::Type{ProteinVariant}, s::AbstractString)::Option{ProteinVariant}
    y = get(_STR_PROTEINVARIANT, s, nothing)
    y === nothing ? none : some(y)
end

struct Protein
    variant::ProteinVariant
    orfs::Vector{UnitRange{UInt16}}
end

struct SegmentData
    gi::Option{String} # missing in human samples
    host::Species
    fluspecies::FluSpecies
    segment::Segment
    subtype::SubType
    year::Int
    name::String
    proteins::Vector{Protein}
    seq::RefValue{Option{LongDNASeq}}
end

function Base.parse(::Type{SubType}, s::Union{String, SubString})::Option{SubType}
    isempty(s) && return none
    s == "Victoria" && return some(SubType(0x00, 0x00, 0x01))
    s == "Yamagata" && return some(SubType(0x00, 0x00, 0x02))
    m = match(r"^H(\d+)N(\d+)$", s)
    m === nothing && error("UNKNOWN SUBTYPE: $s")
    some(SubType(parse(UInt8, m[1]), parse(UInt8, m[2]), 0x00))
end

function parse_cleaned_genomeset(io::IO)::Dict{String, SegmentData}
    result = Dict{String, SegmentData}()
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        data = @unwrap_or parse(SegmentData, line) continue
        gi = unwrap(data.gi)
        @assert !haskey(result, gi) "GB identifier $(gi) not unique"
        result[gi] = data
    end
    result
end

function Base.parse(::Type{SegmentData}, line::Union{String, SubString{String}})::Option{SegmentData}
    fields = split(line, '\t')
    subtype = @? parse(SubType, fields[4])
    year = @? parse_year(fields[6])
    host = Species(fields[2])
    fluspecies = @? parse_fluspecies(fields[8])
    gi = String(fields[1])
    name = String(fields[8])
    segment = @? parse_from_integer(Segment, fields[3])
    some(SegmentData(some(gi), host, fluspecies, segment, subtype, year, name, Protein[], Ref(none(LongDNASeq))))
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

function add_sequences!(segment_data::Dict{String, SegmentData}, io::IO)
    n_updated = 0
    for record in FASTA.Reader(io)
        # headers begins with e.g. gi|59292|gb|X53029|Influenza A virus
        gb_accession = split(FASTA.identifier(record)::String, '|')[4]
        haskey(segment_data, gb_accession) || continue
        data = segment_data[gb_accession]
        @assert is_error(data.seq[])
        data.seq[] = some(FASTA.sequence(LongDNASeq, record))
        n_updated += 1
    end
    println("Updated $n_updated/$(length(segment_data)) records")
    segment_data
end

function parse_inf_aa(io::IO)::Dict{String, ProteinVariant}
    result = Dict{String, ProteinVariant}()
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        fields = split(line, '\t')
        accession = first(fields)
        variant = @unwrap_or parse(ProteinVariant, fields[3]) continue
        @assert !haskey(result, accession) "Duplicate key $accession"
        result[accession] = variant
    end
    result
end


function add_orfs!(segment_data::Dict{String, SegmentData}, accession_protein_map::Dict{String, ProteinVariant}, io::IO)
    n_updates = 0
    proteinbuffer = Protein[]
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        fields = split(line, '\t')
        gb_accession = first(fields)
        haskey(segment_data, gb_accession) || continue
        @assert isodd(length(fields))
        data = segment_data[gb_accession]
        n_updates += 1
        is_bad = false
        empty!(proteinbuffer)
        for (protein_accession, orf_field) in zip(@view(fields[2:2:end]), @view(fields[3:2:end]))

            # For some reason, some of these protein accessions are not actually present
            # in the influenza_aa.dat file.
            if !haskey(accession_protein_map, protein_accession)
                is_bad = true
                break
            end
            variant = accession_protein_map[protein_accession]
            orfs = parse_orf_field(orf_field)
            if is_error(orfs)
                is_bad = true
                break
            end
            push!(proteinbuffer, Protein(variant, unwrap(orfs)))
        end
        if !is_bad
            append!(data.proteins, proteinbuffer)
        end
    end
    println("Updated $n_updates/$(length(segment_data)) records")
    segment_data
end
        
function parse_orf_field(s::Union{String, SubString{String}})::Option{Vector{UnitRange{UInt16}}}
    # If it looks like gb|AB266090:<411->632, the ORF is not present
    # in the reference, and we skip it
    if occursin('>', s) || occursin('<', s)
        return none
    end

    # Multiple ORFs in a gene makes them enclosed in brackets
    # e.g. (gb|AB212651:26-52, 741-1007)
    s1 = strip(s, ['(', ')'])
    p = findfirst(isequal(':'), s1)
    p === nothing && return none
    s2 = @view s1[p+1:ncodeunits(s1)]
    orf_strings = split(s2, ", ")

    orfstrings = split(s2, ", ")
    result = map(parse_range, orfstrings)

    # Return none if it's not divisible by three and thus cant be ORF
    iszero(sum(length, result) % 3) || return none
    some(result)
end

function parse_range(s::AbstractString)::UnitRange{UInt16}
    p2 = findfirst(isequal('-'), s)
    # If it's not a range, it must be a single number
    if p2 === nothing
        n = parse(UInt16, s)
        n:n
    else
        start = parse(UInt16, @view s[1:p2-1])
        stop = parse(UInt16, s[p2+1:ncodeunits(s)])
        start:stop
    end
end

function filter_segment_data!(segment_data::Dict{String, SegmentData})
    len = length(segment_data)
    n_human = count(v -> v.host == human, values(segment_data))

    # Must have a nucleotide sequence
    filter!(segment_data) do (name, data)
        !is_error(data.seq[])
    end
    println("Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # Must be a minimum number of ORFs depending on segment an inf A/B
    filter!(segment_data) do (name, data)
        contains_minimum_proteins(data)
    end
    println("Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # Must not be more than 4 ambiguous bases
    filter!(segment_data) do (name, data)
        count(isambiguous, unwrap(data.seq[])) < 5
    end
    println("Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # Filter for outliers of lengths
    @warn "Not filtering influenza B for segment lengths"
    filter!(segment_data) do (name, data)
        has_acceptable_seq_length(data)
    end
    println("Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # All ORFs can be translated and has stop exactly at end
    filter!(segment_data) do (name, data)
        is_all_translatable(data)
    end
    println("Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    @assert n_human == count(v -> v.host == human, values(segment_data))
    
    return segment_data
end

function contains_minimum_proteins(data::SegmentData)::Bool
    n_proteins = length(data.proteins)
    if data.segment in (PB2, PB1, HA, NP)
        return !iszero(n_proteins)
    end
    if data.segment == NA
        return n_proteins ≥ (data.fluspecies == InfluenzaA ? 1 : 2)
    end
    if data.segment in (MP, NS)
        return n_proteins ≥ 2
    end
    if data.segment == PA
        return n_proteins ≥ (data.fluspecies == InfluenzaA ? 2 : 1)
    end
    @assert false "Unreachable"
end

# This has been empirially determined by looking at the distributions
# in order to find outliers
const ACCEPTABLE_SEGMENT_LENGTHS = Dict(
    NP => 1490:1570,
    HA => 1680:1780,
    MP => 980:1030,
    PB1 => 2270:2345,
    PA => 2150:2240,
    NA => 1350:1470,  # very wide distribution??
    NS => 820:895,
    PB2 => 2280:2345
)

function has_acceptable_seq_length(data::SegmentData)::Bool
    # The sizes for influenza B are a little different.
    data.fluspecies == InfluenzaB && return true
    length(unwrap(data.seq[])) in ACCEPTABLE_SEGMENT_LENGTHS[data.segment]
end

function is_all_translatable(data::SegmentData)::Bool
    seq = unwrap(data.seq[])
    aa_sequence = LongAminoAcidSeq()
    nt_sequence = LongDNASeq()
    for protein in data.proteins
        join!(nt_sequence, (@view(seq[orf]) for orf in protein.orfs))

        # Must have a length divisible by 3
        iszero(length(nt_sequence) % 3) || return false
        BioSequences.translate!(aa_sequence, nt_sequence)
        stop_pos = findfirst(AA_Term, aa_sequence)

        # Must have a stop exactly at the end, nowhere else
        stop_pos == lastindex(aa_sequence) || return false
    end
    return true
end

function add_human_sequences!(segment_data::Dict{String, SegmentData}, io::IO, humanpath::String)
    # First, get a dict of header => seq
    seq_by_header = Dict{String, LongDNASeq}()
    open(humanpath) do io
        reader = FASTA.Reader(GzipDecompressorStream(io))
        for record in reader
            seq_by_header[FASTA.header(record)] = FASTA.sequence(LongDNASeq, record)
        end
    end

    chunks = partition_chunks(io)
    for chunk in chunks
        data = parse_annot_segment_data(chunk)
        @assert !haskey(segment_data, data.name)
        data.seq[] = some(seq_by_header[data.name])
        segment_data[data.name] = data
    end

    segment_data
end

# Please have a look at the files in annotation.tbl to understand this
function partition_chunks(io::IO)::Vector{Vector{String}}
    lines = collect(eachline(io))
    chunks = Vector{Vector{String}}()
    chunk = String[]
    for line in lines
        if startswith(line, "=========")
            push!(chunks, copy(chunk))
            empty!(chunk)
        else
            push!(chunk, line)
        end
    end
    push!(chunks, chunk)
end

function parse_annot_segment_data(chunk::Vector{String})::SegmentData
    @assert startswith(first(chunk), ">Feature ")
    header = first(chunk)[10:end]
    headerfields = split(header, '|')
    species = if startswith(first(headerfields), "A/")
        InfluenzaA
    elseif startswith(first(headerfields), "B/")
        InfluenzaB
    else
        error(first(headerfields))
    end
    subtype = unwrap(parse(SubType, headerfields[2]))
    year = parse(Int, last(split(first(headerfields), '/')))
    segment = expect(parse_from_name(Segment, headerfields[3]), "Error on chunk $header")
    proteins = Protein[]
    orfs = UnitRange{UInt16}[]
    reading_cds = false
    for line in @view chunk[2:end]
        fields = split(line)
        if isnumeric(line[1]) && length(fields) > 2 && fields[3] == "CDS"
            reading_cds = true
        end
        if isnumeric(line[1]) && reading_cds
            push!(orfs, parse(UInt16, fields[1]) : parse(UInt16, fields[2]))
        end
        if reading_cds && fields[1] == "gene"
            variant = expect(parse(ProteinVariant, fields[2]), "Error on chunk $header")
            reading_cds = false
            push!(proteins, Protein(variant, copy(orfs)))
            empty!(orfs)
        end
    end
    SegmentData(none(String), human, species, segment, subtype, year, header, proteins, Ref(none(LongDNASeq)))
end

"Obtain a set of accepted gi for avian and swine"
function cd_hit_deduplicate(segment_data::Dict{String, SegmentData}
)::Dict{Species, Dict{Segment, Set{String}}}
    @warn "Not deduplicating human sequences"

    subdir = joinpath("results", "cdhit")
    isdir(subdir) || mkdir(subdir)

    byspecies = Dict(sp => SegmentData[] for sp in instances(Species))
    for data in values(segment_data)
        push!(byspecies[data.host], data)
    end

    result = Dict{Species, Dict{Segment, Set{String}}}()
    for species in (swine, avian)
        speciesdir = joinpath(subdir, string(species))
        isdir(speciesdir) || mkdir(speciesdir)
        println("Deduplicating $(string(species))...")
        result[species] = cd_hit_deduplicate(byspecies[species], speciesdir)
    end

    # We simply add all human sequences here - no deduplication
    human_dict = Dict(segment => Set{String}() for segment in instances(Segment))
    result[human] = human_dict
    for data in byspecies[human]
        push!(human_dict[data.segment], data.name)
    end

    result
end

function cd_hit_deduplicate(data::Vector{SegmentData}, dirname::String
)::Dict{Segment, Set{String}}
    # Collect FASTAS by segment
    bysegment = Dict(s => FASTA.Record[] for s in instances(Segment))
    for segmentdata in data
        record = FASTA.Record(unwrap(segmentdata.gi), unwrap(segmentdata.seq[]))
        push!(bysegment[segmentdata.segment], record)
    end

    # Write FASTA paths
    for (segment, seqs) in bysegment
        fasta_path = joinpath(dirname, string(segment) * ".fna")
        open(FASTA.Writer, fasta_path) do writer
            for seq in seqs
                write(writer, seq)
            end
        end
    end

    # Deduplicate
    fasta_paths = [joinpath(dirname, string(segment) * ".fna") for segment in instances(Segment)]
    run_cd_hit(fasta_paths)

    # Load in the deduplicated ones
    result = Dict(s => Set{String}() for s in instances(Segment))
    for segment in instances(Segment)
        result[segment] = open(joinpath(dirname, string(segment) * ".fna.cdhit")) do file
            eachline(file) |>
            Filter(x -> startswith(x, '>')) ⨟
            Map(x -> x[2:end]) |>
            Set
        end
    end

    result
end

function run_cd_hit(paths::Vector{String})
    Threads.@threads for path in paths
        command = `bin/cd-hit-est -i $path -o $path.cdhit -aS 0.9 -c 0.95`
        pipe = pipeline(command, stdout="$path.log")
        run(pipe)
    end
end

function serialize_segments(segment_data::Dict{String, SegmentData},
    deduplicated::Dict{Species, Dict{Segment, Set{String}}}
)
    serial_dir = joinpath("results/deduplicated")
    isdir(serial_dir) || mkdir(serial_dir)

    for (species, segment_dict) in deduplicated
        species_dir = joinpath(serial_dir, string(species))
        isdir(species_dir) || mkdir(species_dir)

        # Serialize the ORFs itself
        for (segment, accessions) in segment_dict
            jls_path = joinpath(species_dir, string(segment) * ".jls")
            serialize_orfs(segment_data, accessions, jls_path) # TODO
        end

        # Write the deduplicated FASTA files
        for (segment, accessions) in segment_dict
            fasta_path = joinpath(species_dir, string(segment) * ".fna")
            open(FASTA.Writer, fasta_path) do writer
                for accession in accessions
                    data = segment_data[accession]
                    header = data.host == human ? data.name : unwrap(data.gi)
                    record = FASTA.Record(header, unwrap(data.seq[]))
                    write(writer, record)
                end
            end
        end
    end
end

function serialize_orfs(segment_data::Dict{String, SegmentData},
    accessions::Set{String}, path::String)
    #                 ProteinVar  Vector of orfs
    ProteinType = Tuple{UInt8, Vector{UnitRange{UInt16}}}
    #                accession  vector of proteins
    SegmentType = Tuple{String, Vector{ProteinType}}
    result = Vector{SegmentType}()
    
    for accession in accessions
        data = segment_data[accession]
        vec = Vector{ProteinType}()
        segment_repr = (accession, vec)
        push!(result, segment_repr)
        for protein in data.proteins
            var_uint8 = reinterpret(UInt8, protein.variant)
            push!(vec, (var_uint8, protein.orfs))
        end
    end
    serialize(path, result)
end


isinteractive() || main()
