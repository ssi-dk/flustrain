struct ORFData
    variant::Protein
    seq::LongDNASeq
    aaseq::LongAminoAcidSeq
    orfs::Vector{UnitRange{UInt16}}
    identity::Option{Float64}
    errors::Vector{ErrorMessage}
    # indel errors are special because if there are 50 of them, we don't display them all
    indel_errors::Vector{ErrorMessage}
end

function ORFData(protein::ProteinORF, aln::PairwiseAlignment{LongDNASeq, LongDNASeq},
    ref::Reference
)
    orfseq, orfs, messages, indel_messages = gather_aln(protein, aln)
    aaseq = BioSequences.translate(orfseq)
    ref_aas = (i for (i,n) in zip(ref.seq, protein.mask) if n)
    # Last 3 nts are the stop codon
    refaa = BioSequences.translate(LongDNASeq(collect(ref_aas)[1:end-3]))
    aaaln = pairalign(GlobalAlignment(), aaseq, refaa, AA_ALN_MODEL).aln
    identity = get_identity(aaaln)
    ORFData(protein.var, orfseq, aaseq, orfs, identity, messages, indel_messages)
end

const ALN_MODEL = AffineGapScoreModel(EDNAFULL, gap_open=-25, gap_extend=-2)
const AA_ALN_MODEL = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-2)

is_stop(x::DNACodon) = (x === mer"TAA") | (x === mer"TAG") | (x === mer"TGA")

"Adds one nt at the end of the codon, moving it. If nt is ambiguous, return AAA"
function push_codon(x::DNACodon, nt::DNA)
    val = @inbounds BioSequences.twobitnucs[reinterpret(UInt8, nt) + 1]
    enc = (reinterpret(UInt64, x) << 2 | val) & UInt64(0x3f)
    reinterpret(DNACodon, ifelse(val === 0xff, zero(UInt), enc))
end

# A function to get the Assembly and a protein in a Reference
function gather_aln(protein::ProteinORF, aln::PairwiseAlignment{LongDNASeq, LongDNASeq})
    nucleotides = DNA[]
    last_ref_pos = findlast(protein.mask)
    codon = mer"AAA" # arbitrary starting codon
    seg_pos = ref_pos = n_deletions = n_insertions = 0
    maybe_expected_stop = none(Int)
    # indel messages are special because we choose to skip if there are too many
    messages, indel_messages = ErrorMessage[], ErrorMessage[]
    orfs = UnitRange{UInt16}[]
    orfstart = nothing
    for (seg_nt, ref_nt) in aln
        seg_pos += (seg_nt !== DNA_Gap)
        ref_pos += (ref_nt !== DNA_Gap)
        is_coding = protein.mask[ref_pos]

        if is_coding
            if (orfstart === nothing) & (seg_nt !== DNA_Gap)
                orfstart = seg_pos
            end
        else
            if orfstart !== nothing
                push!(orfs, UInt16(orfstart):UInt16(seg_pos - 1))
                orfstart = nothing

            end
            
            # If we're in an intron, don't do anything
            if ref_pos < last_ref_pos
                continue
            end
        end

        # Check for deletions
        if seg_nt == DNA_Gap
            n_deletions += 1
        else
            codon = push_codon(codon, seg_nt)
            push!(nucleotides, seg_nt)
            if !iszero(n_deletions)
                if !iszero(n_deletions % 3)
                    push!(indel_messages, ErrorMessage(is_important(protein.var) ? important : trivial,
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
                    push!(indel_messages, ErrorMessage(is_important(protein.var) ? important : trivial,
                        "$(protein.var) frameshift insertion of bases $(seg_pos-n_insertions)-$(seg_pos-1)"))
                else
                    push!(indel_messages, ErrorMessage(trivial,
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
                n_aa = div(length(nucleotides), 3)
                expected_aa = div(sum(protein.mask), 3)
                push!(messages, ErrorMessage(is_important(protein.var) ? important : trivial,
                "$(protein.var) stops early, at pos $seg_pos, after $(n_aa) aa (ref is $(expected_aa) aa)."))
            else
                expected_stop = unwrap(maybe_expected_stop)
                if expected_stop != seg_pos
                    @assert seg_pos > expected_stop
                    push!(messages, ErrorMessage(trivial,
                    "$(protein.var) stops late, at pos $seg_pos, expected pos $expected_stop."))
                end
            end
            break
        end
    end

    if orfstart !== nothing
        push!(orfs, UInt16(orfstart):UInt16(seg_pos))
    end

    dnaseq = LongDNASeq(nucleotides)
    # Is seq length divisible by three?
    remnant = length(dnaseq) % 3
    if !iszero(remnant)
        push!(messages, ErrorMessage(is_important(protein.var) ? important : trivial,
            "$(protein.var) length is $(length(dnaseq)) nt, not divisible by 3, last nt trimmed off"))
    end

    # Does it end with a stop? If so, remove it, else report error
    dnaseq = iszero(remnant) ? dnaseq : dnaseq[1:end-remnant]
    if isempty(dnaseq) || !is_stop(DNACodon(dnaseq[end-2:end]))
        push!(messages, ErrorMessage(is_important(protein.var) ? important : trivial,
            "$(protein.var) does not end with a stop codon."))
    else
        dnaseq = dnaseq[1:end-3]
    end

    return dnaseq, orfs, messages, indel_messages
end