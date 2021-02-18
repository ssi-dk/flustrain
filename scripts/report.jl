# TODO: Check for extra nt at ends, by aligning to a selection of conserved ends

using FASTX
using BioSequences
using BioAlignments
using CodecZlib
using ErrorTypes
using Printf
using Serialization
using Plots

# Number of nucleotides from each end to consider the "ends" of the sequence
# in which we have more lax requirements for depth
const TERMINAL = 25

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

"""
    Segment(::String) -> Segment

Create a representaiton of one of the eight Influenza A genome segments.
"""
@enum Segment::UInt8 PB1 PB2 PA HA NP NA MP NS
const _STR_SEGMENT = Dict(map(i -> string(i)=>i, instances(Segment)))
function Segment(s::AbstractString)
    smt = get(_STR_SEGMENT, s, nothing)
    smt === nothing && error("Unknown segment: \"$s\"")
    smt
end

@enum ProteinVariant::UInt8 pb2 pb1 pb1fa pa pax ha np na m1 m2 ns1 nep
function Base.print(io::IO, var::ProteinVariant)
    str = if var == pb1fa
        "PB1-FA"
    elseif var == pax
        "PA-X"
    else
        uppercase(string(Symbol(var)))
    end
    print(io, str)
end

const _CRITICAL = Bool[1,1,0,1,0,1,1,1,1,1,1,1]
is_critical(x::ProteinVariant) = @inbounds _CRITICAL[reinterpret(UInt8, x) + 0x01]

const ALN_MODEL = AffineGapScoreModel(EDNAFULL, gap_open=-25, gap_extend=-2)

# Data retrieved from the .mat.gz file
struct Depths
    depths::Vector{UInt32}
end

# Data retrieved from the reference files
struct Protein
    var::ProteinVariant
    # a bitvector over its reference seq with 1s for coding nucleotides, incl last stop
    mask::BitVector
end

struct Reference
    seq::LongDNASeq
    proteins::Vector{Protein}
end

# Data retrieved from the .fsa files
struct Assembly
    insignificant::BitVector
    seq::LongDNASeq
    accession::String
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
    proteins::Vector{Protein}
end

struct ErrorMessage
    critical::Bool
    msg::String
end

function merge_data_sources(assemblies::Dict{Segment, Option{Assembly}},
    depths::Dict{Segment, Option{Vector{UInt32}}},
    references::Dict{Tuple{Segment, String}, Reference},
    asm_identities::Dict{Segment, Option{Float64}}
)::Vector{Tuple{Segment, Option{SegmentData}}}
    
    result = Vector{Tuple{Segment, Option{SegmentData}}}()
    for (segment, maybe_assembly) in assemblies
        if is_none(maybe_assembly)
            push!(result, (segment, none))
            continue
        end
        depth = expect(depths[segment],
            "Segment $segment of basename $basename found in assembly, but not in matrix"
        )
        assembly = unwrap(maybe_assembly)
        reference = references[(segment, assembly.accession)]
        asm_id = unwrap(asm_identities[segment])
        (seq, ref) = (assembly.seq, reference.seq)
        aln = pairalign(OverlapAlignment(), seq, ref, ALN_MODEL).aln
        smt = SegmentData(depth, asm_id, assembly.insignificant, seq, ref, aln, reference.proteins)
        push!(result, (segment, Thing(smt)))
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
function load_depths(matpath::AbstractString)::Dict{Segment, Option{Vector{UInt32}}}
    open(matpath) do io
        result = Dict{Segment, Option{Vector{UInt32}}}()
        segment = nothing
        depths = UInt32[]
        fields = Vector{SubString{String}}(undef, 7)
        linedepths = Vector{UInt32}(undef, 6)
        for line in Iterators.map(strip, eachline(GzipDecompressorStream(io)))
            if isempty(line)
                if !isempty(depths)
                    haskey(result, segment) && error("Segment $segment present twice in $matpath")
                    result[segment::Segment] = Thing(depths)
                    depths = UInt32[]
                end
                continue
            end
            if startswith(line, '#')
                # Headers look like "HEADER_HA"
                p = findlast(isequal('_'), line)
                p === nothing && error("Found header \"$(line)\", expected pattern \"HEADER_SEGMENT\"")
                segment = Segment(line[p+1:end])
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
        for segment in instances(Segment)
            haskey(result, segment) || (result[segment] = none)
        end
        result
    end
end

# Here, we load template identities from kma2.res. If a template has <99.5% identity,
# it means the first assembly did not converge.
function load_res_file(resfilename::AbstractString)::Dict{Segment, Option{Float64}}
    open(resfilename) do io
        fields = Vector{SubString{String}}(undef, 11)
        lines = eachline(io)
        header, _ = iterate(lines)::NTuple{2, Any}
        result = Dict(s => ErrorTypes.None{Float64}() for s in instances(Segment))
        @assert header == "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value"
        for fields in (lines |> imap(strip) |> ifilter(!isempty) |> imap(x -> split!(fields, x, UInt8('\t'))))
            segment = Segment(strip(last(rsplit(first(fields), '_', limit=2))))
            is_none(result[segment]) || error("Segment $(string(segment)) present twice in file $resfilename")
            result[segment] = Thing(parse(Float64, fields[5]) / 100)
        end
        result
    end
end

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
                accession = FASTA.identifier(record)
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
                    reference = Reference(seq, Protein[])
                    result[(segment, accession)] = reference
                    for (var_uint8, orf_tuple) in v
                        mask = falses(length(seq))
                        for orf in orf_tuple
                            @view(mask[orf]) .= true
                        end
                        push!(reference.proteins, Protein(ProteinVariant(var_uint8), mask))
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

"Given an aln b/w ref and whole segment, return (seqnucs, refnucs, pos_in_seg)"
function gather_aln(protein::Protein, aln::PairwiseAlignment{LongDNASeq, LongDNASeq})
    seg_nucs, ref_nucs, seg_positions = DNA[], DNA[], UInt16[]
    seg_pos = ref_pos = num_seg_nt = zero(UInt16)
    last_ref_pos = findlast(protein.mask)
    codon = mer"AAA"
    expected_stop = ErrorTypes.None{Int}()
    for (seg_nt, ref_nt) in aln
        seg_pos += (seg_nt !== DNA_Gap)
        ref_pos += (ref_nt !== DNA_Gap)
        is_coding = protein.mask[ref_pos]

        # If we're in an intron, don't do anything
        !is_coding & (ref_pos < last_ref_pos) && continue

        if seg_nt !== DNA_Gap
            num_seg_nt += UInt16(1)
            codon = push_codon(codon, seg_nt)
        end

        if ref_pos == last_ref_pos
            expected_stop = Thing(seg_pos % Int)
        end

        push!(seg_nucs, seg_nt)
        push!(ref_nucs, ref_nt)
        push!(seg_positions, seg_pos)

        # Only stop if we find a stop codon NOT in an intron
        iszero(num_seg_nt % 3) & is_stop(codon) && break
    end
    return LongDNASeq(seg_nucs), LongDNASeq(ref_nucs), seg_positions, expected_stop
end

function protein_errors(protein::Protein, aln::PairwiseAlignment{LongDNASeq, LongDNASeq}
)::Tuple{Vector{ErrorMessage}, Vector{ErrorMessage}}
    seg_nucs, ref_nucs, seg_positions, expected_stop = gather_aln(protein, aln)
    errors, indel_messages = ErrorMessage[], ErrorMessage[]
    n_del = n_ins = 0
    for (seg_nt, ref_nt, pos) in zip(seg_nucs, ref_nucs, seg_positions)
        if seg_nt == DNA_Gap
            n_del += 1
        elseif !iszero(n_del)
            if !iszero(n_del % 3)
                push!(indel_messages, ErrorMessage(is_critical(protein.var),
                    "$(protein.var) frameshift deletion of $n_del bases b/w bases $(pos-1)-$(pos)"))
            else
                push!(indel_messages, ErrorMessage(false,
                    "$(protein.var) deletion of $n_del bases b/w bases $(pos-1)-$(pos)"))
            end
            n_del = 0
        end

        if ref_nt == DNA_Gap
            n_ins += 1
        elseif !iszero(n_ins)
            if !iszero(n_ins % 3)
                push!(indel_messages, ErrorMessage(is_critical(protein.var),
                    "$(protein.var) frameshift insertion of bases $(pos-n_ins)-$(pos-1)"))
            else
                push!(indel_messages, ErrorMessage(false,
                    "$(protein.var)  insertion of bases $(pos-n_ins)-$(pos-1)"))
            end
            n_ins = 0
        end
    end
    seq = ungap(seg_nucs)
    
    # Is seq length divisible by three?
    remnant = length(seq) % 3
    if !iszero(remnant)
        push!(errors, ErrorMessage(is_critical(protein.var),
            "$(protein.var) length is $(length(seq)) nt, not divisible by 3"))
    end

    # Does it not end with a stop codon?
    aaseq = BioSequences.translate(iszero(remnant) ? seq : seq[1:end-remnant])
    stoppos = findfirst(AA_Term, aaseq)
    if stoppos === nothing
        push!(errors, ErrorMessage(is_critical(protein.var),
            "$(protein.var) does not have a stop codon."))
    end

    # Does it have a stop codon too early or too late?
    if is_none(expected_stop)
        push!(errors, ErrorMessage(is_critical(protein.var),
        "$(protein.var) has an early stop at pos $(last(seg_positions))."))
    else
        exp_stop = unwrap(expected_stop)
        if exp_stop != last(seg_positions)
            @assert last(seg_positions) > exp_stop
            push!(errors, ErrorMessage(false,
            "$(protein.var) stops late, at pos $(last(seg_positions)), expected pos $(exp_stop)."))
        end
    end
    
    sort!(errors, rev=true, by=x -> x.critical)
    sort!(indel_messages, rev=true, by=x -> x.critical)
    return (errors, indel_messages)
end

# This function gets the identity between a segment and its reference
function get_identity(aln::PairwiseAlignment{LongDNASeq, LongDNASeq})
    seq, ref = fill(DNA_Gap, length(aln)), fill(DNA_Gap, length(aln))
    for (i, (seqnt, refnt)) in enumerate(aln)
        seq[i] = seqnt
        ref[i] = refnt
    end
    start = max(findfirst(!isgap, seq), findfirst(!isgap, ref))
    stop = min(findlast(!isgap, seq), findlast(!isgap, ref))
    id = count(zip(view(seq, start:stop), view(ref, start:stop))) do pair
        first(pair) === last(pair)
    end / length(start:stop)
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

function load_assembly(assemblypath::AbstractString)::Dict{Segment, Option{Assembly}}
    open(FASTA.Reader, assemblypath) do reader
        result = Dict{Segment, Option{Assembly}}()
        map(reader) do record
            id = FASTA.identifier(record)::String
            v = rsplit(id, '_', limit=2)
            length(v) == 2 || error("Found header \"$(id)\", expected pattern \"HEADER_SEGMENT\"")
            accession = first(v)
            segment = Segment(last(v))
            seq = FASTA.sequence(LongDNASeq, record)
            @assert length(record.sequence) == length(seq)
            insignificant = BitVector([in(i, UInt8('a'):UInt8('z')) for i in @view record.data[record.sequence]])
            haskey(result, segment) && error("Segment $segment present twice in $assemblypath")
            result[segment] = Thing(Assembly(insignificant, seq, accession))
        end
        for segment in instances(Segment)
            haskey(result, segment) || (result[segment] = none)
        end
        result
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
    if is_none(maybe_data)
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
    idstr = "identity $(@sprintf "%.3f" identity)"
    push!(lines, beginning * " $depthstr $covstr $idstr")

    # Warning if coverage or identity is too low
    coverage < 0.9 && push!(lines, "\t\t ERROR: Coverage is less than 90%")
    identity < 0.9 && push!(lines, "\t\t ERROR: Identity is less than 90%")

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

    for protein in data.proteins
        (errors, indel_messages) = protein_errors(protein, data.aln)

        if length(indel_messages) > 4
            push!(lines, "\t\t" * (is_critical(protein.var) ? "ERROR:" : "      ") *
            " Numerous indels in $(protein.var)")
        else
            for msg in indel_messages
                push_msg(lines, msg, 2)
            end
        end
        for msg in errors
            push_msg(lines, msg, 2)
        end
    end
    lines
end

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
        for segment_dict in assemblies
            for (segment, maybe_assembly) in segment_dict
                assembly = @unwrap_or maybe_assembly continue
                push!(accessions, (segment, assembly.accession))
            end
        end
        load_references(accessions, refdir)
    end

    # Now merge it
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

    # Create consensus files
    for (basename, data_vec) in segment_data
        open(FASTA.Writer, joinpath(outdir, basename, "consensus.fna")) do writer
            for (segment, maybe_data) in data_vec
                data = @unwrap_or maybe_data continue
                write(writer, FASTA.Record(basename, segment, data.seq, data.insignificant))
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
