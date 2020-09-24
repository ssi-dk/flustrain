import tools
from collections import Counter

def getskip(str, i):
    n = []
    for i in range(i, len(str)):
        c = str[i]
        if not c.isdigit():
            break
        n.append(c)

    return len(n), int(''.join(n))

def get_consensus(pos, str, ref):
    kinds = Counter()
    insertions = Counter()
    skip = 0
    for (i, c) in enumerate(str):
        if skip > 0:
            skip -= 1
            continue

        if c in ("*", "#"):
            kinds[''] += 1
        elif c == '^':
            skip = 1
        elif c in ('>', '<', '$'):
            pass
        elif c == '-':
            n, s = getskip(str, i+1)
            skip = n + s
        elif c == '+':
            n, s = getskip(str, i+1)
            skip = n + s
            insertions[str[i+n+1:i+n+1+s]] += 1
        elif c in ('.', ','):
            kinds[ref] += 1
        else:
            assert c.upper() in {'A', 'C', 'G', 'T'}
            kinds[c.upper()] += 1

    ncalls = sum(kinds.values())
    base, basecount = kinds.most_common()[0]

    if ncalls < 10:
        base = ref
    elif insertions:
        ins, N = insertions.most_common()[0]
        if N > ncalls / 2:
            base = ins + base

    return (base, ncalls)

def construct_consensus(pairs, filename):
    "pairs of (base, ncalls)"
    buffer = list()
    for (i, (b, n)) in enumerate(pairs):
        if i > 20 and i < len(pairs) - 20 and n < 10:
            raise ValueError("Depth drops to {} at pos {} in {}".format(n, i, filename))
        buffer.append(b)

    return ''.join(buffer)

def parse_consensus(filename, filehandle, subject):
    pieces = list()
    reflen = 0
    for line in filehandle:
        fields = line.split()
        if fields[0] != subject:
            continue
        reflen += 1
        ref = fields[2]
        pos = int(fields[1])
        pieces.append(get_consensus(pos, fields[4], ref))

    seq = construct_consensus(pieces, filename)

    if reflen < 50:
        raise ValueError("Reference gene in {} less than 50 bp. Error?".format(filename))

    if len(seq) < 0.9 * reflen:
        raise ValueError("Sequence {} less than 90% of reference".format(filename))

    return reflen, seq

def translate_and_verify(subject, sequence, reflen):
    if len(sequence) % 3 != 0:
        sequence = sequence[:-(len(sequence) % 3)]

    aa = tools.translate(sequence)

    # Find first stop, complain if it's too early and truncate.
    stoppos = aa.find("*")
    if stoppos > -1:
        if stoppos < 0.9 * len(sequence):
            msg = ("Stoppos at position {}/{} in {} "
                  "more than 10% of sequence truncated.").format(stoppos+1, len(aa), snakemake.input[0])
            raise ValueError(msg)
        else:
            aa = aa[:stoppos]

    return aa

with open(snakemake.input[1]) as file:
    subject = next(file).strip()

with open(snakemake.input[0]) as file:
    reflen, nt = parse_consensus(str(snakemake.input[0]), file, subject)

aa = translate_and_verify(subject, nt, reflen)

with open(snakemake.output[0], "w") as file:
    print(">{}\n{}".format(snakemake.wildcards.basename, nt), file=file)

with open(snakemake.output[1], "w") as file:
    print(">{}\n{}".format(snakemake.wildcards.basename, aa), file=file)
