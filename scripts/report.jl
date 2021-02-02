using FASTX
using BioSequences
using BioAlignments
using CodecZlib
using ErrorTypes

const TERMINAL = 25

"""
    Segment(::String) -> Segment

Create a representaiton of one of the eight Influenza A genome segments.
"""
@enum Segment::UInt8 begin
    PB1
    PB2
    PA
    HA
    NP
    NA
    MP
    NS
end

const _STR_SEGMENT = Dict(map(i -> string(i)=>i, instances(Segment)))
Segment(s::AbstractString) = _STR_SEGMENT[s]

# Splicing together the nucleotides indexed by the frames gives the mature
# ORF which can be translated to the full-length protein.
# Non-critical proteins may be truncated in viable virus.
struct Protein
    name::String
    segment::Segment
    frames::Vector{UnitRange{Int}}
    critical::Bool
end

const PROTEINS = Protein[
    Protein("HA", HA, [17:1717], true),
    Protein("M1", MP, [14:769], true),
    Protein("M2", MP, [14:39, 728:992], true),
    Protein("NA", NA, [9:1418], true),
    Protein("NP", NP, [34:1527], true),
    Protein("NS1", NS, [15:665], true),
    Protein("NEP", NS, [15:44, 517:849], true),
    Protein("PA", PA, [13:2160], true),
    Protein("PAX", PA, [13:582, 584:769], false),
    Protein("PB1", PB1, [13:2283], true),
    Protein("PB1-F2", PB1, [221:376], false),
    Protein("PB2", PB2, [16:2292], true)
]

const PROTEINS_OF_SEGMENT = Dict(s => Protein[] for s in instances(Segment))
for protein in PROTEINS
    push!(PROTEINS_OF_SEGMENT[protein.segment], protein)
end

#const REFERENCES = open(FASTA.Reader, joinpath(dirname(abspath(@__FILE__)), "orf_ref.fna")) do reader
const REFERENCES = open(FASTA.Reader, joinpath(dirname(abspath(".")), "orf_ref.fna")) do reader
    res = Dict{Segment, LongDNASeq}()
    for record in reader
        segment = Segment(FASTA.identifier(record))
        sequence = FASTA.sequence(LongDNASeq, record)
        res[segment] = sequence
    end
    res
end

const MODEL_1 = AffineGapScoreModel(EDNAFULL, gap_open=-13, gap_extend=-2)

function expected_length(p::Protein)
    nt = sum(length, p.frames)
    @assert iszero(rem(nt, 3))
    div(nt, 3)
end

function joinseqs(v::AbstractVector{T}) where T <: LongSequence
    ln = sum(length, v)
    seq = T(ln)
    offset = 0
    for s in v
        seq[offset+1:offset+length(s)] = s
        offset += length(s)
    end
    seq
end

struct Error
    msg::String
end

function errormessage(x::Result{T, Error}) where T
    is_error(x) || error("Attempted error message call on non-error")
    x.data._1.msg
end

"Given the orf of `ref`, extract the equivalent seq in `query`"
function extract_orf(ref::LongDNASeq, query::LongDNASeq)::Result{LongDNASeq, Error}
    # OverlapAlignment is critical here, the others fucks it up
    aln = pairalign(OverlapAlignment(), query, ref, MODEL_1).aln

    # Alignment must cover all of the reference
    if aln.a.aln.firstref != 1 || aln.a.aln.lastref != lastindex(ref)
        return Err(Error("Sequence does not align to ORF"))
    end

    # Get the area from the first ref align to last ref align
    seqaln, refaln = Vector{DNA}(undef, length(aln)), Vector{DNA}(undef, length(aln))
    for (i, (q, r)) in enumerate(aln)
        seqaln[i] = q
        refaln[i] = r
    end

    leading_gaps = findfirst(!isgap, refaln) - 1
    trailing_gaps = findfirst(!isgap, reverse!(refaln)) - 1

    result = ungap(LongDNASeq(@view(seqaln[1 + leading_gaps : end - trailing_gaps])))

    if length(result) < 0.9 * length(ref)
        msg = "Extracted orf too short, expected $(length(ref)), got $(length(result))"
        return Err(Error(msg))
    end

    return Ok(result)
end

"Join and extract a sequence of one or more ORFs"
function extract_seq(refseq::LongDNASeq, orfs::Vector{<:UnitRange}, query::LongDNASeq
)::Result{LongSequence, Error}
    fragments = LongDNASeq[]
    for orf in orfs
        push!(fragments, @?(extract_orf(refseq[orf], query)))
    end
    Ok(joinseqs(fragments))
end

function get_depth(path::AbstractString)::Option{Tuple{Vector{UInt32}, LongDNASeq}}
    reference = DNA[]
    depths = UInt32[]

    isfile(path) || return none
    open(path) do io
        lines = eachline(GzipDecompressorStream(io))
        iterate(lines) === nothing && return none # skip header
        stripped = (rstrip(line) for line in lines)
        for fields in (split(strip) for strip in stripped if !isempty(strip))
            nucleotide = convert(DNA, first(first(fields)))
            nucleotide == DNA_Gap && continue

            # A, C, G, T and also -, but not N
            depth = sum(x -> parse(UInt32, x), @view fields[2:5]) + parse(UInt32, fields[7])
            push!(depths, depth)
            push!(reference, nucleotide)
        end
        Thing((depths, LongDNASeq(reference)))
    end
end

# Depth naturally declines towards the ends of the template. This is not
# necessarily a sign of bad sequencing. So I include the first and last 25 bp,
# and check for depth < 25.
count_low_depths(x::Vector{<:Integer}) = count(<(25), @view x[1+TERMINAL:end-TERMINAL])

"Returns a (is_accepted, error_string)"
function protein_errors(query::LongDNASeq, protein::Protein)::Tuple{Bool, Option{String}}
    # Get the sequence corresponding to the protein, if possible.
    extracted = extract_seq(REFERENCES[protein.segment], protein.frames, query)
    is_error(extracted) && return (false, Thing("Error: $(errormessage(extracted))"))
    extracted_seq = unwrap(extracted)

    if !iszero(rem(length(extracted_seq), 3))
        return (false, Thing("Error: Length of ORF $(protein.name) not divisible by three"))
    end

    aas = translate(extracted_seq)
    termpos = findfirst(AA_Term, aas)
    len = termpos === nothing ? length(aas) : termpos - 1

    expected = expected_length(protein)
    if protein.critical && len < 0.9 * expected
        return (false, Thing("Error: Length of ORF $(protein.name) less than 90% of ref"))
    elseif len != expected
        factor = round(len / expected, digits=2)
        return (true, Thing("Warning: Length of ORF $(protein.name) $(factor) % of ref"))
    end
    (true, none)
end

"Loads sequence and number of lowercase letters"
function load_consensus(string::AbstractString)::Option{Tuple{LongDNASeq, UInt}}
    isfile(string) || return none
    reader = FASTA.Reader(open(string))
    it = iterate(reader)
    it === nothing && return none
    rec, _ = it
    seq = FASTA.sequence(LongDNASeq, rec)
    n_insig = count(byte -> islowercase(Char(byte)), @view rec.data[rec.sequence]) % UInt
    Thing((seq, n_insig))
end

####################################
"Return is_accepted, errors/warnings"
function report_text(conspath::AbstractString, matpath::AbstractString, segment::Segment
)::Tuple{Bool, Vector{String}}

    # Get depths and ref from mat, return stirng if none
    dep = get_depth(matpath)
    if is_none(dep)
        return (false, ["Error: Empty or missing matrix file."])
    end
    depths, refseq = unwrap(dep)

    # Get consens from conspath, return string if none
    con = load_consensus(conspath)
    if is_none(con)
        return (false, ["Error: Empty or missing consensus file"])
    end
    consensus, n_insignificant = unwrap(con)

    errors = String[]
    # Check for insignificant bases
    if !iszero(n_insignificant)
        push!(errors, ["Warning: $n_insignificant bases insignificantly basecalled"])
    end

    # Check for low depths
    lowdepth = count_low_depths(depths)
    if !iszero(lowdepth)
        push!(errors, ["Warning: $lowdepth bases with depth < 25 (disregarding terminals)"])
    end

    # Check for ambiguous bases
    amb = count(isambiguous, consensus)
    if !iszero(amb)
        push!(errors, ["Warning: $amb ambiguous bases in consensus sequence"])
    end

    # For each protein, check it.
    for protein in PROTEINS_OF_SEGMENT[segment]
        prot_error = protein_errors(consensus, protein)
        if !is_none(prot_error)
            push!(errors, unwrap(prot_error))
        end
    end

    return (true, errors)
end

# Check all files are unique and all segments

####################################

struct ConsensusState
    segment::String
    ok::Bool
    mean_depth::Float64
    coverage::Float64
    n_low_depth::Int
    n_insignificant::Int
    n_ns::Int
    orf_lens::Tuple{Int, Int} # observed / ref
end

function ConsensusState(segment::AbstractString, conspath::AbstractString,
    matpath::AbstractString
)::ConsensusState
    # If file is empty
    first_it = open(iterate, FASTA.Reader, conspath)
    EMPTY_STATE = ConsensusState(segment, false, 0.0, 0.0, 0, 0, 0)
    isnothing(first_it) && return EMPTY_STATE

    consensus::FASTA.Record = first(first_it)
    isempty(consensus.sequence) && return EMPTY_STATE

    depths, reference = get_depth(matpath)
    n_low_depth = count_low_depths(depths)
    n_ns = count_ns(consensus)
    n_insignificant = count_insignificant_bases(consensus)
    mean_depth = sum(depths) / length(depths)
    coverage = count(!iszero, depths) / length(depths)

    ok = iszero(n_low_depth) && iszero(n_ns) && iszero(n_insignificant)

    return ConsensusState(segment, ok, mean_depth, coverage, n_low_depth, n_insignificant, n_ns)
end

function report_text(state::ConsensusState)::String
    dep = round(state.mean_depth, digits=2)
    cov = round(state.coverage, digits=2)
    lines = ["$(state.segment):\tDepth: $dep\tCoverage: $cov"]
    # If the state is OK, we just return that.
    state.ok && return first(lines)

    if !iszero(state.n_low_depth)
        push!(lines, "\t$(state.n_low_depth) bases (excluding terminal $(TERMINAL)bp) has depth < 10")
    end
    if !iszero(state.n_ns)
        push!(lines, "\t$(state.n_ns) bases called as \"N\"")
    end
    if !iszero(state.n_insignificant)
        push!(lines, "\t$(state.n_insignificant) basecalls without statistical significance")
    end
    join(lines, '\n')
end

function main(conspaths::Vector{<:AbstractString}, matpath::Vector{<:AbstractString},
    accpath::AbstractString, reportpath::AbstractString
)
    # Step one: Calculate the report info for each segment
    states = ConsensusState[]
    for (conspath, matpath) in zip(conspaths, matpaths)
        segment = first(splitext(basename(conspath)))
        push!(states, ConsensusState(segment, conspath, matpath))
    end

    # Then print the accepted ones:
    open(accpath, "w") do accfile
        for state in states
            state.ok && println(accfile, state.segment)
        end
    end

    # Then print report text
    open(reportpath, "w") do reportfile
        for state in states
            println(reportfile, report_text(state), '\n')
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 4
        error("Usage: julia report.jl conspaths matpaths accpath reportpath")
    end
    conspaths = split(ARGS[1], ',')
    matpaths = split(ARGS[2], ',')
    length(conspaths) == length(matpaths) || error("must have same number of cons and matpaths")
    main(conspaths, matpaths, ARGS[3], ARGS[4])
end
