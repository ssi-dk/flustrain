using FASTX
using BioSequences
using BioAlignments
using CodecZlib
using ErrorTypes

const TERMINAL = 25

struct Protein
    name::String
    segment::String
    frames::Vector{UnitRange{Int}}
    critical::Bool
end

const PROTEINS = Protein[
    Protein("HA", "HA", [17:1717], true),
    Protein("M1", "MP", [14:769], true),
    Protein("M2", "MP", [14:39, 728:992], true),
    Protein("NA", "NA", [9:1418], true),
    Protein("NP", "NP", [34:1527], true),
    Protein("NS1", "NS", [15:665], true),
    Protein("NEP", "NS", [15:44, 517:849], true),
    Protein("PA", "PA", [13:2160], true),
    Protein("PAX", "PA", [13:582, 584:769], false),
    Protein("PB1", "PB1", [13:2283], true),
    Protein("PB1-F2", "PB1", [221:376], false),
    Protein("PB2", "PB2", [16:2292], true)
]

#const REFERENCES = open(FASTA.Reader, joinpath(dirname(abspath(@__FILE__)), "orf_ref.fna")) do reader
const REFERENCES = open(FASTA.Reader, joinpath(dirname(abspath(".")), "orf_ref.fna")) do reader
    res = Dict{String, LongDNASeq}()
    for record in reader
        segment = FASTA.identifier(record)
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

"Given the orf of `ref`, extract the equivalent seq in `query`"
function extract_orf(ref::LongDNASeq, query::LongDNASeq)
    # OverlapAlignment is critical here, the others fucks it up
    aln = pairalign(OverlapAlignment(), query, ref, MODEL_1).aln

    # Alignment must cover all of the reference
    @assert aln.a.aln.firstref == 1
    @assert aln.a.aln.lastref == lastindex(ref)

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
        error("Extracted orf too short, expected $(length(ref)), got $(length(result))")
    end

    return result
end

"Join and extract a sequence of one or more ORFs"
function extract_seq(refseq::LongDNASeq, orfs::Vector{<:UnitRange}, query::LongDNASeq)
    joinseqs([extract_orf(refseq[orf], query) for orf in orfs])
end

function validate_protein(query::LongDNASeq, protein::Protein)
    accepted = true
    errors = String[]
    extracted = extract_seq(REFERENCES[protein.segment], protein.frames, query)

    if !iszero(rem(length(extracted), 3))
        accepted = false
        push!(errors, "ERROR: Length of ORF $(protein.name) not divisible by three")
        return accepted, errors
    end

    aas = translate(extracted)
    termpos = findfirst(x -> x == AA_Term, collect(aas))
    len = termpos === nothing ? length(aas) : termpos - 1

    expected = expected_length(protein)
    if protein.critical && len < 0.9 * expected
        accepted = false
        push!(errors, "ERROR: Length of ORF $(protein.name) less than 90% of ref")
    elseif len != expected
        factor = round(len / expected, digits=2)
        push!(warnings, "WARNING: Length of ORF $(protein.name) $(factor) % of ref")
    end

    return accepted, errors
end

function get_depth(path::AbstractString)
    reference = DNA[]
    depths = UInt32[]

    open(path) do io
        lines = eachline(GzipDecompressorStream(io))
        iterate(lines) # skip header
        stripped = (rstrip(line) for line in lines)
        for fields in (split(strip) for strip in stripped if !isempty(strip))
            nucleotide = convert(DNA, first(first(fields)))
            nucleotide == DNA_Gap && continue

            # A, C, G, T and also -, but not N
            depth = sum(x -> parse(UInt32, x), @view fields[2:5]) + parse(UInt32, fields[7])
            push!(depths, depth)
            push!(reference, nucleotide)
        end
        (depths, LongDNASeq(reference))
    end
end

# Depth naturally declines towards the ends of the template. This is not
# necessarily a sign of bad sequencing. So I include the first and last 25 bp,
# and check for depth < 25.
count_low_depths(x::Vector{<:Integer}) = count(<(25), @view x[1+TERMINAL:end-TERMINAL])

function count_insignificant_bases(rec::FASTA.Record)
    count(byte -> islowercase(Char(byte)), @view rec.data[rec.sequence])
end

count_ns(rec::FASTA.Record) = count(@view rec.data[rec.sequence]) do byte
    (byte === UInt8('n')) | (byte === UInt8('N'))
end

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
