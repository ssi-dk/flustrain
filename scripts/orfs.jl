using FASTX
using BioSequences
using Setfield

struct ORF
    part::UnitRange{UInt32}

    function ORF(part::UnitRange{UInt32})
        isempty(part) && error("Cannot instantiate empty ORF")
        iszero(length(part) % 3) || error("Length of part not divisible by three: $part")
        new(part)
    end
end

codon(::DNAAlphabet) = DNACodon
codon(::DNAAlphabet) = RNACodon
function isstop(x::Union{DNACodon, RNACodon})
    (x === reinterpret(typeof(x), mer"TAG")) |
    (x === reinterpret(typeof(x), mer"TAA")) |
    (x === reinterpret(typeof(x), mer"TGA"))
end
isstart(x::Union{DNACodon, RNACodon}) = x === reinterpret(typeof(x), mer"ATG")

# TODO: This only works with met starts
# Also, does not call overlapping ORFs
function forward_orfs(seq::BioSequence{<:NucleicAcidAlphabet{2}})
    orfs = ORF[]
    isempty(seq) && return orfs
    starts = (0, 0, 0)
    local mer
    local start
    for outer mer in each(codon(Alphabet(seq)), seq)
        frame = mer.position % 3 + 1
        start = starts[frame]
        if isstop(mer.fw)
            iszero(start) && continue
            orf = ORF(UInt32(start):UInt32(mer.position) - UInt32(1))
            push!(orfs, orf)
            starts = (@set starts[frame] = 0)
        elseif isstart(mer.fw)
            newstart = iszero(start) ? mer.position : min(start, mer.position)
            starts = (@set starts[frame] = newstart)
        end
    end
    # Last ORF, which does not end with STOP
    if !isstop(mer.fw) && !iszero(start)
        orf = ORF(UInt32(start):UInt32(mer.position) - UInt32(1))
        push!(orfs, orf)
    end
    orfs
end

    for frame in 1:3
        stop = length(seq) - (length(seq) - frame + 1) % 3
        aa_seq = translate(seq[frame:stop])



#=


GENE_ORFS = [
    ("HA", "HA", [17:1720]),
    ("M1", "MP", [14:772]),
    ("M2", "MP", [14:39, 728:995]),
    ("NA", "NA", [9:1421]), # 1:1419 for N5
    ("NP", "NP", [34:1530]),
    ("NS1", "NS", [15:668]),
    ("NEP", "NS", [15:44, 517:852]),
    ("PA", "PA", [13:2163]),
    ("PAX", "PA", [13:582, 584:772]),
    ("PB1", "PB1", [13:2286]),
    ("PB1-F2", "PB1", [221:379]),
    ("PB2", "PB2", [16:2295])
]



=#