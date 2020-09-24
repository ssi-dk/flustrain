import os
import gzip
import tools
from collections import defaultdict

TERMINAL = 10

def get_depth(filename):
    "Return depth, N low coverage pos (not incl 10 terminal), reference"
    # Get by-base depth
    depths = list()
    reference = list()
    with gzip.open(filename, "rt") as file:
        next(file) # skip header
        for fields in filter(None, map(lambda line: line.rstrip().split(), file)):
            # Skip bases that are gap in reference (i.e. insertions)
            if fields[0] == '-':
                continue

            # A,C,G or T
            depths.append(sum(map(int, fields[1:5])))
            reference.append(fields[0].upper())

    lowcov = sum(1 for i in depths[TERMINAL:-TERMINAL] if i < 10)
    refseq = tools.FastaEntry("ref", bytearray(''.join(reference).encode()))
    return sum(depths) / len(depths), lowcov, refseq

report_file = open(snakemake.output[0], "w")
accepted_file = open(snakemake.output[1], "w")
consensuses = dict()
subtypes = dict()
accepted = dict()

kinds = [snakemake.input.con, snakemake.input.mat, snakemake.input.sub]
for (conpath, matpath, subpath) in zip(*[sorted(i) for i in kinds]):
    gene, _ = os.path.splitext(os.path.basename(conpath))
    with open(subpath) as file:
        subtype = next(file).strip()

    with open(conpath, "rb") as file:
        consensus = next(tools.byte_iterfasta(file))

    consensuses[gene] = consensus
    subtypes[gene] = subtype
    accepted[gene] = True

    print(f"{gene}:\t{subtypes[gene]} ", file=report_file, end="")

    if len(consensus) == 0:
        continue

    depth, low_coverage, reference = get_depth(matpath)
    print(f"depth: {depth:>9.2f}", file=report_file)
    # Check low coverage
    if low_coverage > 0:
        print(f"\t{low_coverage}/{len(reference)-2*TERMINAL} bases (excluding terminals) has depth < 10", file=report_file)
        accepted[gene] = False

    # Check length
    if len(consensus) != len(reference):
        _, orfc = tools.find_orf(consensus)
        _, orfr = tools.find_orf(reference)
        if abs(len(orfc) - len(orfr)) > 12:
            print(f"\tReference ORF len: {len(orfr)} consensus ORF len: {len(orfc)}", file=report_file)
            accepted[gene] = False

    middleNs = consensus.sequence[TERMINAL:-TERMINAL].upper().count(b'N')
    if middleNs > 0:
        print(f"\t{middleNs} Ns excluding terminal bp", file=report_file)
        accepted[gene] = False

    Ns = consensus.sequence.upper().count(b'N')
    if Ns > 9:
        print(f"\t{Ns} Ns total", file=report_file)
        accepted[gene] = False

    # Check lowercase (not accurate bases)
    lowercases = (ord('a'), ord('c'), ord('g'), ord('t'))
    inaccurate = sum(1 for i in consensus.sequence[TERMINAL:-TERMINAL] if i in lowercases)
    if inaccurate > 0:
        print(f"\tInaccurate bases (excluding terminals): {inaccurate}/{len(consensus.sequence) - 2*TERMINAL}", file=report_file)
        accepted[gene] = False

# Check that all subtypes is the same:
subtype_set = set(subtypes.values())
if len(set(subtype_set)) > 1:
    bysubtype = defaultdict(list)
    for (gene, subtype) in subtypes.items():
        bysubtype[subtype].append(gene)
    print("\nDIFFERING SUBTYPES:", file=report_file)
    for (subtype, genes) in bysubtype.items():
        print(f"\t{subtype}: {','.join(genes)}", file=report_file)

    for gene in subtypes:
        accepted[gene] = False

print("", file=report_file)

for (gene, acc) in sorted(accepted.items()):
    if acc:
        print(gene, file=accepted_file)

accepted_file.close()
report_file.close()
