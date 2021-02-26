# TODO: Check for extra nt at ends, by aligning to a selection of conserved ends

module t

using FASTX
using BioSequences
using BioAlignments
using CodecZlib
using ErrorTypes
using Influenza
using Printf
using Serialization
using Plots

# Number of nucleotides from each end to consider the "ends" of the sequence
# in which we have more lax requirements for depth
const TERMINAL = 25

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

const SegmentTuple{T} = NTuple{length(instances(Segment)), T}

const _CRITICAL = Tuple(Bool[1,1,0,0,1,0,1,1,1,1,1,1,1,0,1,1,1])
@assert length(_CRITICAL) == length(instances(Protein))
is_critical(x::Protein) = @inbounds _CRITICAL[reinterpret(UInt8, x) + 0x01]

const ALN_MODEL = AffineGapScoreModel(EDNAFULL, gap_open=-25, gap_extend=-2)
const AA_ALN_MODEL = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-2)

# Data retrieved from the .mat.gz file
struct Depths
    depths::Vector{UInt32}
end

# Data retrieved from the reference files
struct ProteinORF
    var::Protein
    # a bitvector over its reference seq with 1s for coding nucleotides, incl last stop
    mask::BitVector
end

struct Reference
    seq::LongDNASeq
    proteins::Vector{ProteinORF}
end

# Data retrieved from the .fsa files
struct Assembly
    insignificant::BitVector
    seq::LongDNASeq
    accession::String
end

struct ErrorMessage
    critical::Bool
    msg::String
end

struct ORFData
    variant::Protein
    seq::LongDNASeq
    aaseq::LongAminoAcidSeq
    identity::Option{Float64}
    errors::Vector{ErrorMessage}
    # indel errors are special because if there are 50 of them, we don't display them all
    indel_errors::Vector{ErrorMessage}
end

# Synthesized from all of the above. Its construction should verify that they are
# consistent with each other (which is distinct from them being internally consistent)
struct SegmentData
    depths::Vector{UInt32}
    assembly_identity::Float64
    insignificant::BitVector
    seq::LongDNASeq
    refseq::LongDNASeq
    aln::PairwiseAlignment{LongDNASeq, LongDNASeq}
    orfdata::Vector{ORFData}
end

"Create SegmentData from all data sources below"
function merge_data_sources(assemblies::SegmentTuple{Option{Assembly}},
    depths::SegmentTuple{Option{Vector{UInt32}}},
    references::Dict{Tuple{Segment, String}, Reference},
    asm_identities::SegmentTuple{Option{Float64}}
)::Vector{Tuple{Segment, Option{SegmentData}}}
    
    result = Vector{Tuple{Segment, Option{SegmentData}}}()
    for (segment_index, segment) in enumerate(instances(Segment))
        if is_error(assemblies[segment_index])
            push!(result, (segment, none))
            continue
        end
        depth = expect(depths[segment_index],
            "Segment $segment of basename $basename found in assembly, but not in matrix"
        )
        assembly = unwrap(assemblies[segment_index])
        reference = references[(segment, assembly.accession)]
        asm_id = unwrap(asm_identities[segment_index])
        (seq, ref) = (assembly.seq, reference.seq)
        aln = pairalign(OverlapAlignment(), seq, ref, ALN_MODEL).aln
        orfdata = map(reference.proteins) do protein
            ORFData(protein, aln, reference)
        end
        smt = SegmentData(depth, asm_id, assembly.insignificant, seq, ref, aln, orfdata)
        push!(result, (segment, some(smt)))
    end
    sort!(result, by=first)
end

function Base.join(v::AbstractVector{T}) where T <: LongSequence
    ln = sum(length, v)
    seq = T(ln)
    offset = 0
    @inbounds for s in v
        seq[offset+1:offset+length(s)] = s
        offset += length(s)
    end
    seq
end

function split!(v::Vector{SubString{String}}, s::Union{String, SubString{String}}, sep::UInt8)
    n = 0
    start = 1
    @inbounds for i in 1:ncodeunits(s)
        if codeunit(s, i) == sep
            n += 1
            n >= length(v) && throw(BoundsError(v, n+1))
            substr = SubString(s, start, i-1)
            v[n] = substr
            start = i + 1
        end
    end
    @inbounds v[n+1] = SubString(s, start, ncodeunits(s))
    v
end

"Given a path to a .mat.gz file, return a dict of an optional vec of depths"
function load_depths(matpath::AbstractString)::SegmentTuple{Option{Vector{UInt32}}}
    open(matpath) do io
        result = fill(none(Vector{UInt32}), length(instances(Segment)))
        segment = nothing
        depths = UInt32[]
        fields = Vector{SubString{String}}(undef, 7)
        linedepths = Vector{UInt32}(undef, 6)
        for line in eachline(GzipDecompressorStream(io)) |> imap(string)
            if isempty(line)
                if !isempty(depths)
                    segment_index = reinterpret(UInt8, segment::Segment) + 0x01
                    is_error(result[segment_index]) || error("Segment $segment present twice in $matpath")
                    result[segment_index] = some(depths)
                    depths = UInt32[]
                end
                continue
            end
            if startswith(line, '#')
                # Headers look like "HEADER_HA"
                p = findlast(isequal('_'), line)
                p === nothing && error("Found header \"$(line)\", expected pattern \"HEADER_SEGMENT\"")
                segment = unwrap(parse(Segment, line[p+1:end]))
                continue
            end
            split!(fields, line, UInt8('\t'))
            @inbounds for i in 1:6
                linedepths[i] = parse(UInt32, fields[i+1])
            end
            # Skip lines with majority as deletion
            argmax(linedepths) == 6 && continue
            total_depth = sum(linedepths)
            # We don't include N in depths!
            depth = total_depth - linedepths[5]
            push!(depths, depth)
        end
        SegmentTuple(result)
    end
end

# Here, we load template identities from kma2.res. If a template has <99.5% identity,
# it means the first assembly did not converge.
function load_res_file(resfilename::AbstractString)::SegmentTuple{Option{Float64}}
    open(resfilename) do io
        fields = Vector{SubString{String}}(undef, 11)
        lines = eachline(io)
        header, _ = iterate(lines)::NTuple{2, Any}
        result = fill(none(Float64), length(instances(Segment)))
        @assert header == "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value"
        for fields in (lines |> imap(strip) |> ifilter(!isempty) |> imap(x -> split!(fields, x, UInt8('\t'))))
            segment = unwrap(parse(Segment, strip(last(rsplit(first(fields), '_', limit=2)))))
            segment_index = reinterpret(UInt8, segment) + 0x01
            is_error(result[segment_index]) || error("Segment $(string(segment)) present twice in file $resfilename")
            result[segment_index] = some(parse(Float64, fields[5]) / 100)
        end
        SegmentTuple(result)
    end
end

"Given a set of (segment, accession), fetches References for each. A Reference is the sequence
and its Proteins"
function load_references(references::Set{Tuple{Segment, String}}, refdir::AbstractString
)::Dict{Tuple{Segment, String}, Reference}
    # Group segments
    bysegment = Dict{Segment, Set{String}}()
    for (segment, accession) in references
        push!(get!(Set{String}, bysegment, segment), accession)
    end

    record = FASTA.Record()
    result = Dict{Tuple{Segment, String}, Reference}()
    for (segment, accessions) in bysegment
        seqpath = joinpath(refdir, "$segment.fna")
        orfpath = joinpath(refdir, "$segment.orfs.jls")

        # Get the sequence
        seqof = Dict{String, LongDNASeq}()
        open(FASTA.Reader, seqpath) do reader
            while !eof(reader)
                read!(reader, record)
                accession = FASTA.identifier(record)::String
                if in(accession, accessions)
                    seqof[accession] = FASTA.sequence(LongDNASeq, record)
                end
            end
        end

        # Get the proteins
        added_accessions = Set{String}()
        open(orfpath) do io
            for (accession, v) in deserialize(io)
                if in(accession, accessions)
                    if !haskey(seqof, accession)
                        error("Accession $accession missing from $seqpath")
                    end
                    seq = seqof[accession]
                    reference = Reference(seq, ProteinORF[])
                    result[(segment, accession)] = reference
                    for (var_uint8, orf_tuple) in v
                        mask = falses(length(seq))
                        for orf in orf_tuple
                            @view(mask[orf]) .= true
                        end
                        push!(reference.proteins, ProteinORF(Protein(var_uint8), mask))
                    end
                    push!(added_accessions, accession)
                end
            end
        end
        
        # Check we added all accessions
        missing_acc = setdiff(accessions, added_accessions)
        if !isempty(missing_acc)
            error("Accession $(first(missing_acc)) missing from $orfpath")
        end
    end
    result
end

function push_msg(lines::Vector{<:AbstractString}, x::ErrorMessage, indent::Integer)
    push!(lines, "$('\t'^indent)$(x.critical ? "ERROR: " : "       ")$(x.msg)")
end

is_stop(x::DNACodon) = (x === mer"TAA") | (x === mer"TAG") | (x === mer"TGA")

"Adds one nt at the end of the codon, moving it. If nt is ambiguous, return AAA"
function push_codon(x::DNACodon, nt::DNA)
    val = @inbounds BioSequences.twobitnucs[reinterpret(UInt8, nt) + 1]
    enc = (reinterpret(UInt64, x) << 2 | val) & UInt64(0x3f)
    reinterpret(DNACodon, ifelse(val === 0xff, zero(UInt), enc))
end

# A function to get the Assembly and the Reference
function gather_aln(protein::ProteinORF, aln::PairwiseAlignment{LongDNASeq, LongDNASeq})
    nucleotides = DNA[]
    last_ref_pos = findlast(protein.mask)
    codon = mer"AAA" # arbitrary starting codon
    seg_pos = ref_pos = n_deletions = n_insertions = 0
    maybe_expected_stop = none(Int)
    # indel messages are special because we choose to skip if there are too many
    messages, indel_messages = ErrorMessage[], ErrorMessage[]
    for (seg_nt, ref_nt) in aln
        seg_pos += (seg_nt !== DNA_Gap)
        ref_pos += (ref_nt !== DNA_Gap)
        is_coding = protein.mask[ref_pos]

        # If we're in an intron, don't do anything
        !is_coding & (ref_pos < last_ref_pos) && continue

        # Check for deletions
        if seg_nt == DNA_Gap
            n_deletions += 1
        else
            codon = push_codon(codon, seg_nt)
            push!(nucleotides, seg_nt)
            if !iszero(n_deletions)
                if !iszero(n_deletions % 3)
                    push!(indel_messages, ErrorMessage(is_critical(protein.var),
                        "$(protein.var) frameshift deletion of $n_deletions bases b/w bases $(seg_pos-1)-$(seg_pos)"))
                else
                    push!(indel_messages, ErrorMessage(false,
                        "$(protein.var) deletion of $n_deletions bases b/w bases $(seg_pos-1)-$(seg_pos)"))
                end
                n_deletions = 0
            end
        end

        # Check for insertions
        if ref_nt == DNA_Gap
            n_insertions += 1
        else
            if !iszero(n_insertions)
                if !iszero(n_insertions % 3)
                    push!(indel_messages, ErrorMessage(is_critical(protein.var),
                        "$(protein.var) frameshift insertion of bases $(seg_pos-n_insertions)-$(seg_pos-1)"))
                else
                    push!(indel_messages, ErrorMessage(false,
                    "$(protein.var)  insertion of bases $(seg_pos-n_insertions)-$(seg_pos-1)"))
                end
                n_insertions = 0
            end
        end

        if ref_pos == last_ref_pos
            maybe_expected_stop = some(seg_pos % Int)
        end

        # Only stop if we find a stop codon NOT in an intron
        if is_stop(codon) && iszero(length(nucleotides) % 3)
            if is_error(maybe_expected_stop)
                push!(messages, ErrorMessage(is_critical(protein.var),
                "$(protein.var) has an early stop at segment pos $seg_pos, after $(length(nucleotides)) nt."))
            else
                expected_stop = unwrap(maybe_expected_stop)
                if expected_stop != seg_pos
                    @assert seg_pos > expected_stop
                    push!(messages, ErrorMessage(false,
                    "$(protein.var) stops late, at pos $seg_pos, expected pos $expected_stop."))
                end
            end
            break
        end
    end

    dnaseq = LongDNASeq(nucleotides)
    # Is seq length divisible by three?
    remnant = length(dnaseq) % 3
    if !iszero(remnant)
        push!(messages, ErrorMessage(is_critical(protein.var),
            "$(protein.var) length is $(length(dnaseq)) nt, not divisible by 3, last nt trimmed off"))
    end

    # Does it end with a stop? If so, remove it, else report error
    dnaseq = iszero(remnant) ? dnaseq : dnaseq[1:end-remnant]
    if isempty(dnaseq) || !is_stop(DNACodon(dnaseq[end-2:end]))
        push!(messages, ErrorMessage(is_critical(protein.var),
            "$(protein.var) does not end with a stop codon."))
    else
        dnaseq = dnaseq[1:end-3]
    end

    return dnaseq, messages, indel_messages
end

function ORFData(protein::ProteinORF, aln::PairwiseAlignment{LongDNASeq, LongDNASeq},
    ref::Reference)
    orfseq, messages, indel_messages = gather_aln(protein, aln)
    aaseq = BioSequences.translate(orfseq)
    ref_aas = (i for (i,n) in zip(ref.seq, protein.mask) if n)
    # Last 3 nts are the stop codon
    refaa = BioSequences.translate(LongDNASeq(collect(ref_aas)[1:end-3]))
    aaaln = pairalign(GlobalAlignment(), aaseq, refaa, AA_ALN_MODEL).aln
    identity = get_identity(aaaln)
    ORFData(protein.var, orfseq, aaseq, identity, messages, indel_messages)
end

# This function gets the identity between a segment and its reference
function get_identity(aln::PairwiseAlignment{T, T})::Option{Float64} where {T <: BioSequence}
    iszero(length(aln)) && return none

    gapval = gap(eltype(T))
    seq, ref = fill(gapval, length(aln)), fill(gapval, length(aln))
    for (i, (seqnt, refnt)) in enumerate(aln)
        seq[i] = seqnt
        ref[i] = refnt
    end

    # This can sometimes happen if e.g. there is no called sequence
    # in the supposed amino acid sequence.
    (all(isgap, seq) || all(isgap, ref)) && return none

    start = max(findfirst(!isgap, seq), findfirst(!isgap, ref))
    stop = min(findlast(!isgap, seq), findlast(!isgap, ref))
    id = count(zip(view(seq, start:stop), view(ref, start:stop))) do pair
        first(pair) === last(pair)
    end / length(start:stop)

    return some(id)
end

# We only have this to reconstruct the record - including the lowercase letters.
function FASTA.Record(basename::AbstractString, segment::Segment,
    seq::LongDNASeq, insignificant::BitVector
)
    data = UInt8['>'; codeunits(basename); '_'; codeunits(string(segment)); '\n']
    seqdata = take!((io = IOBuffer(); print(io, seq); io))
    @assert length(seqdata) == length(insignificant)
    for (i, is_insignificant) in enumerate(insignificant)
        is_insignificant && (seqdata[i] += 0x20)
    end
    FASTA.Record(append!(data, seqdata))
end

function load_assembly(assemblypath::AbstractString)::SegmentTuple{Option{Assembly}}
    open(FASTA.Reader, assemblypath) do reader
        result = fill(none(Assembly), length(instances(Segment)))
        foreach(reader) do record
            id = FASTA.identifier(record)::String
            v = rsplit(id, '_', limit=2)
            length(v) == 2 || error("Found header \"$(id)\", expected pattern \"HEADER_SEGMENT\"")
            accession = first(v)
            segment = unwrap(parse(Segment, last(v)))
            segment_index = reinterpret(UInt8, segment) + 0x01
            seq = FASTA.sequence(LongDNASeq, record)
            @assert length(record.sequence) == length(seq)
            insignificant = BitVector([in(i, UInt8('a'):UInt8('z')) for i in @view record.data[record.sequence]])
            is_error(result[segment_index]) || error("Segment $segment present twice in $assemblypath")
            result[segment_index] = some(Assembly(insignificant, seq, accession))
        end
        SegmentTuple(result)
    end
end

function report_lines(basename::AbstractString, data_vec::Vector{Tuple{Segment, Option{SegmentData}}})
    lines = [basename * ':']
    for (segment, maybe_data) in data_vec
        push_report_lines!(lines, segment, maybe_data)
    end
    push!(lines, "")
end

function push_report_lines!(lines::Vector{<:AbstractString}, segment, maybe_data::Option{SegmentData})
    beginning = '\t' * rpad(string(segment) * ":", 4)

    # Check presence of segment
    if is_error(maybe_data)
        push!(lines, beginning * "\n\t\tERROR: Missing segment")
        return nothing
    end
    data = unwrap(maybe_data)

    # Check segment length and return early if way too short
    minlen = min(length(data.seq), length(data.depths))
    if minlen < 2 * TERMINAL + 1
        push!(lines, beginning * "\n\t\tERROR: Sequence or depths length: $minlen")
        return nothing
    end

    # Header: HA: depth 4.32e+03 coverage 1.000 
    mean_depth = sum(UInt, data.depths) / length(data.depths)
    coverage = count(!iszero, data.depths) / length(data.depths)
    identity = get_identity(data.aln)
    depthstr = "depth $(@sprintf "%.2e" mean_depth)"
    covstr = "coverage $(@sprintf "%.3f" coverage)"
    idstr = is_error(identity) ? "N/A" : "identity $(@sprintf "%.3f" unwrap(identity))"
    push!(lines, beginning * " $depthstr $covstr $idstr")

    # Warning if coverage or identity is too low
    coverage < 0.9 && push!(lines, "\t\t ERROR: Coverage is less than 90%")
    is_error(identity) || unwrap(identity) < 0.9 && push!(lines, "\t\t ERROR: Identity is less than 90%")

    # Warning if identity b/w first and second assembly too low (asm not converged)
    if data.assembly_identity < 0.995
        push!(lines, "\t\t       Id. b/w first and second asm low: $(@sprintf "%.4f" data.assembly_identity)")
    end

    # insignificant bases
    n_insignificant = count(data.insignificant)
    if !iszero(n_insignificant)
        push!(lines, "\t\t       $(n_insignificant) bases insignificantly basecalled")
    end

    # Check for low depths
    lowdepth = count(x -> x < 25, @view data.depths[TERMINAL + 1: end - TERMINAL])
    if !iszero(lowdepth)
        push!(lines, "\t\t       $lowdepth central bases with depth < 25")
    end

    # Check for ambiguous bases
    amb = count(isambiguous, data.seq)
    if !iszero(amb)
        push!(lines, "\t\t       $amb ambiguous bases in consensus sequence")
    end

    for orfdata in data.orfdata
        if length(orfdata.indel_errors) > 4
            push!(lines, "\t\t" * (is_critical(orfdata.variant) ? "ERROR:" : "      ") *
            " Numerous indels in $(orfdata.variant)")
        else
            for msg in orfdata.indel_errors
                push_msg(lines, msg, 2)
            end
        end
        for msg in orfdata.errors
            push_msg(lines, msg, 2)
        end
    end
    lines
end

# This function is thread-unsafe, and must be called in serial.
function plot_depths(path::AbstractString, data_vec::Vector{Tuple{Segment, Option{SegmentData}}})
    plt = plot(ylabel="Log10 depths", xticks=nothing, ylim=(-0.1, 5))
    for (segment, maybe_data) in data_vec
        data = @unwrap_or maybe_data continue
        ys = log10.(data.depths)
        xs = range(0.0, stop=1.0, length=length(ys))
        plot!(plt, xs, ys, label=string(segment), legend=:outertopright)
    end
    savefig(plt, path)
end

"Load all input files, consolidating it to SegmentData objects"
function load_all_data(basenames::Vector{String}, alnpath::AbstractString,
    assemblyfilename::AbstractString,  resfilename::AbstractString,
    matfilename::AbstractString, refdir::AbstractString
)::Vector{Tuple{String, Vector{Tuple{Segment, Option{SegmentData}}}}}
    
    # Load depths, assemblies and assembly identities
    depths = Vector{Any}(nothing, length(basenames))
    assemblies = copy(depths)
    identities = copy(depths)

    Threads.@threads for (i, basename) in collect(enumerate(basenames))
        depths[i] = load_depths(joinpath(alnpath, basename, matfilename))
        assemblies[i] = load_assembly(joinpath(alnpath, basename, assemblyfilename))
        identities[i] = load_res_file(joinpath(alnpath, basename, resfilename))
    end

    references::Dict{Tuple{Segment, String}, Reference} = let
        accessions = Set{Tuple{Segment, String}}()
        for segment_vals in assemblies
            for (i, maybe_assembly) in enumerate(segment_vals)
                assembly = @unwrap_or maybe_assembly continue
                push!(accessions, (Segment(i - 1), assembly.accession))
            end
        end
        load_references(accessions, refdir)
    end

    # Now merge the sources of data
    data_vec = Vector{Any}(nothing, length(basenames))
    Threads.@threads for i in eachindex(depths)
        data_vec[i] = merge_data_sources(assemblies[i], depths[i], references, identities[i])
    end
    sort!(collect(zip(basenames, data_vec)), by=first)
end

function main(outdir::AbstractString, reportpath::AbstractString, depthsdir::AbstractString,
    alnpath::AbstractString, assemblyfilename::AbstractString, resfilename::AbstractString,
    matfilename::AbstractString, refdir::AbstractString
)
    basenames = readdir(alnpath)
    segment_data = load_all_data(basenames, alnpath, assemblyfilename, resfilename, matfilename, refdir)

    # Create the report
    lines_vectors = Vector{Vector{String}}(undef, length(basenames))
    Threads.@threads for (i, (basename, data_vec)) in collect(enumerate(segment_data))
        lines_vectors[i] = report_lines(basename, data_vec)
    end

    open(reportpath, "w") do report
        for linevec in lines_vectors
            for line in linevec
                println(report, line)
            end
        end
    end

    # Create consensus files - both DNA and AA
    for (basename, data_vec) in segment_data
        open(FASTA.Writer, joinpath(outdir, basename, "consensus.fna")) do writer
            for (segment, maybe_data) in data_vec
                data = (@unwrap_or maybe_data continue)::SegmentData
                write(writer, FASTA.Record(basename, segment, data.seq, data.insignificant))
            end
        end

        open(FASTA.Writer, joinpath(outdir, basename, "consensus.faa")) do writer
            for (segment, maybe_data) in data_vec
                data = (@unwrap_or maybe_data continue)::SegmentData
                for orfdata in data.orfdata
                    write(writer, FASTA.Record(basename * '_' * string(orfdata.variant), orfdata.aaseq))
                end
            end
        end
    end

    # Create depths plots - plotting in inherently not threadsafe
    for (basename, data_vec) in segment_data
        plot_depths(joinpath(depthsdir, basename * ".pdf"), data_vec)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 8
        error("Usage: julia report.jl outdir reportpath depthsdir alnpath asmname resname matname refdir")
    end
    main(ARGS...)
end

end # module

