"""
    ORFData

An `ORFData` object contains protein-specific information from the alignment
between an influenza segment and a reference segment. Its construction checks for
sutff lifethe presence of valid reading frames, indels compared to the reference.

Construct with `ORFData(::ProteinORF, ::PairwiseAlignment)`.
"""
struct ORFData
    variant::Protein
    seq::LongDNASeq
    aaseq::LongAminoAcidSeq
    orfs::Vector{UnitRange{UInt16}}
    identity::Option{Float64}
    errors::Vector{ErrorMessage}
end

function ORFData(protein::ProteinORF, aln::PairwiseAlignment{LongDNASeq, LongDNASeq},
    ref::Reference
)
    orfseq, orfs, messages = gather_aln(protein, aln)
    aaseq = BioSequences.translate(orfseq)
    ref_aas = (i for (i,n) in zip(ref.seq, protein.mask) if n)
    # Last 3 nts are the stop codon
    refaa = BioSequences.translate(LongDNASeq(collect(ref_aas)[1:end-3]))
    aaaln = pairalign(GlobalAlignment(), aaseq, refaa, DEFAULT_AA_ALN_MODEL).aln
    identity = alignment_identity(aaaln)
    ORFData(protein.var, orfseq, aaseq, orfs, identity, messages)
end

is_stop(x::DNACodon) = (x === mer"TAA") | (x === mer"TAG") | (x === mer"TGA")

"Adds one nt at the end of the codon, moving it. If nt is ambiguous, return AAA"
function push_codon(x::DNACodon, nt::DNA)
    val = @inbounds BioSequences.twobitnucs[reinterpret(UInt8, nt) + 1]
    enc = (reinterpret(UInt64, x) << 2 | val) & UInt64(0x3f)
    reinterpret(DNACodon, ifelse(val === 0xff, zero(UInt), enc))
end

"""Compares a protein and an alignment between a segment containing the protein
and the referece segment that contains that protein.
"""
function gather_aln(protein::ProteinORF, aln::PairwiseAlignment{LongDNASeq, LongDNASeq}
)::Tuple{LongDNASeq, Vector{UnitRange{UInt16}}, Vector{ErrorMessage}}
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
        is_coding = !iszero(ref_pos) && protein.mask[ref_pos]

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
                    severity = if is_important(protein.var) && ref_pos < last_ref_pos - 10
                        important
                    else
                        trivial
                    end
                    push!(indel_messages, ErrorMessage(severity,
                        "$(protein.var) frameshift deletion of $n_deletions bases b/w bases $(seg_pos-1)-$(seg_pos)"))
                else
                    # We accept 7 aas deleted, but not more
                    push!(indel_messages, ErrorMessage(n_deletions > 21 ? important : trivial,
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
                    severity = if is_important(protein.var) && ref_pos < last_ref_pos - 10
                        important
                    else
                        trivial
                    end
                    push!(indel_messages, ErrorMessage(severity,
                        "$(protein.var) frameshift insertion of bases $(seg_pos-n_insertions)-$(seg_pos-1)"))
                else
                    # 12 aa insertions are OK, more must be a nonfunctional mutant
                    push!(indel_messages, ErrorMessage(n_insertions > 36 ? important : trivial,
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

                # We consider early stopping (truncation) an error if it is 15 aa or more shorter
                importance = is_important(protein.var) && n_aa < expected_aa - 14 ? important : trivial
                push!(messages, ErrorMessage(importance,
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

    # Merge indel_msg and messages here
    # If there are too many indel messages, we collapse them
    if length(indel_messages) > 4
        push!(messages, ErrorMessage(is_important(protein.var) ? important : trivial,
            "Numerous (>4) indels in $(protein.var)"))
    else
        append!(messages, indel_messages)
    end

    return dnaseq, orfs, messages
end
