import re
import tools
from collections import defaultdict

class Mutation:
    _delpattern = re.compile("^Δ(\d+)(?:-(\d+))?$")
    _singlepattern = re.compile("^([A-Z])(\d+)([A-Z])$")
    _multipattern = re.compile("^([A-Z](?:/[A-Z])*)(\d+)([A-Z](?:/[A-Z])*)$")
    __slots__ = ['pos', 'ref', 'alt']

    def __init__(self, pos, ref, alt):
        self.pos = int(pos)
        for (n, v) in (("ref", ref), ("alt", alt)):
            if isinstance(v, int):
                v_ = chr(v)
            elif v is None:
                v_ = None
            elif isinstance(v, bytearray):
                v_ = v.decode()
            else:
                v_ = v

            if not (v is None or isinstance(v_, str)):
                raise ValueError(f"Ref and alt must be strings, None, integers or bytearrays, not {type(v_)}")
            self.__setattr__(n, v_)

    @classmethod
    def fromstring(cls, st, reference=None):
        if st.startswith("Δ"):
            m = cls._delpattern.match(st)
            if m is None:
                raise ValueError(f"String {st} cannot be interpreted as mutation")

            start, stop = m.groups()
            start = int(start) - 1
            stop = int(stop)
            if stop is None:
                stop = start
            return cls(start, reference[start:stop], None)

        m = cls._singlepattern.match(st)
        if m is None:
            raise ValueError(f"String {st} cannot be interpreted as mutation")

        ref, pos, alt = m.groups()
        return cls(int(pos)-1, ref, alt)

    @classmethod
    def multistring(cls, st, reference=None):
        m = cls._multipattern.match(st)
        if m is None:
            raise ValueError(f"String {st} cannot be interpreted as (multi)mutation")

        ref, pos, alt = m.groups()
        p = int(pos) - 1
        mutations = list()
        for r in ref.split('/'):
            for a in alt.split('/'):
                mutations.append(cls(p, r, a))

        return mutations

    def __hash__(self):
        return hash(self.ref) ^ hash(self.pos) ^ hash(self.alt)

    def __repr__(self):
        if self.alt is None:
            pend = self.pos + len(self.ref)
            if pend != self.pos:
                return f"Δ{self.pos+1}-{pend}"
            else:
                return f"Δ{self.pos+1}"
        else:
            return f"{self.ref}{self.pos+1}{self.alt}"

    def __eq__(self, other):
        return self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

def remove_dual_gaps(A, B):
    A_, B_ = list(), list()
    for (i, j) in zip(A, B):
        if (i, j) != (45, 45):
            A_.append(i)
            B_.append(j)

    return bytearray(A_).upper(), bytearray(B_).upper()

def get_mutations(reference, sequence):
    if len(reference) != len(sequence):
        raise ValueError("Length of reference and sequence must be identical")

    reference, sequence = remove_dual_gaps(reference, sequence)

    p = -1
    gapsize = 0
    mutations = list()
    for (i, (ref, alt)) in enumerate(zip(reference, sequence)):
        if ref != '-':
            p += 1

        if alt == '-':
            gapsize += 1

        elif gapsize > 0:
            mutations.append(Mutation(p - gapsize, reference[i-gapsize:i], None))
            gapsize = 0

        if ref == '-':
            gapsize -= 1

        elif gapsize < 0:
            mutations.append(Mutation(p-1, reference[i+gapsize-1], sequence[i+gapsize:i]))
            gapsize = 0

        if alt != ref and alt != '-' and ref != '-':
            mutations.append(Mutation(p, ref, alt))

    if gapsize > 0:
        mutations.append(Mutation(p - gapsize, reference[i-gapsize:i], None))
    elif gapsize < 0:
        mutations.append(Mutation(p-1, reference[i+gapsize-1], sequence[i+gapsize:i]))

    return set(mutations)

def get_keys(string, ref):
    if '/' in string and '+' in string:
        raise ValueError("This part not implemented yet!")
    if '/' in string:
        return [frozenset([m]) for m in Mutation.multistring(string, ref)]
    elif '+' in string:
        s = set()
        for mstring in string.split('+'):
            s.add(Mutation.fromstring(mstring, ref))
        return [frozenset(s)]
    else:
        return [frozenset([Mutation.fromstring(string, ref)])]

# {frozenset(Mut1, Mut2): ["NI", "RI" ... ]}
def extract_mutations(reference, table):
    # Split table
    result = dict()
    for row in table.splitlines():
        fields = row.split()
        keys = get_keys(fields[0], reference)
        for key in keys:
            result[key] = fields[1:]

    return result

_H1N1_TABLE = """I117R   NI  RI  NI  NI
E119A   RI  RI  RI  RI
E119D   RI  HRI HRI HRI
E119G   NI  HRI HRI HRI
E119V   RI  HRI RI  NI
Q136K/Q NI  RI  NI  NI
Q136K   NI  HRI HRI RI
Q136R   NI  HRI HRI RI
D151D/E NI  RI  RI  NI
D151N/D RI  RI  NI  NI
R152K   RI  NI  NI  NI
D199E   RI  NI  NI  NI
D199G   RI  NI  NI  NI
I223K   RI  NI  NI  NI
I223R   RI  RI  NI  NI
I223V   NI  NI  NI  NI
I223T   RI  NI  NI  NI
S247N   NI  NI  NI  NI
S247G   RI  NI  NI  NI
S247R   RI  RI  HRI HRI
H275Y   HRI NI  HRI NI
R293K   RI  NI  NI  NI
N295S   HRI NI  RI  NI
D199N+H275Y HRI NI  HRI NI
G147R+H275Y HRI NI  HRI NI
E119A+H275Y HRI RI  NI  NI
E119D+H275Y HRI HRI HRI HRI
E119G+H275Y HRI HRI HRI HRI
I223K+H275Y HRI RI  HRI RI
I223R+H275Y HRI RI  HRI RI
I223V+H275Y HRI NI  HRI NI
S247N+H275Y HRI NI  HRI NI
Q313K+I427T RI  RI  NI  NI"""

_H3N2_TABLE = """E119D   NI  RI  NI  NI
E119I   HRI RI  NI  NI
E119V   HRI NI  NI  NI
Q136K   NI  HRI NI  NI
D151A   NI  RI  RI  NI
D151E   RI  NI  NI  NI
D151G   NI  HRI NI  NI
D151V/D NI  HRI NI  NI
I222L   NI  NI  NI  NI
R224K   HRI RI  NI  NI
Δ245-248   HRI RI  NI  NI
Δ247-250   HRI RI  NI  NI
K249E   RI  NI  NI  NI
E276D   RI  HRI NI  NI
R292K   HRI HRI HRI NI
N294S   HRI NI  NI  NI
N329K    RI  RI  NI  NI
S331R    RI  RI  NI  NI
R371K   RI  RI  NI  NI
Q391K   RI  RI  RI  NI
E119V+T148I HRI HRI HRI HRI
E119V+I222L HRI RI  NI  NI
E119V+I222V HRI NI  NI  NI
I222T+S331R RI  NI  NI  NI"""

# Returns [ (refmutset, resistance) tuples ]
def get_resistances(reference, consensus, subtype):
    if subtype == "H1N1":
        table = extract_mutations(reference, _H1N1_TABLE)
    elif subtype == "H3N2":
        table = extract_mutations(reference, _H3N2_TABLE)
    else:
        raise ValueError(f"No mutations known for subtype {subtype}")

    mutations = get_mutations(reference, consensus)
    result = list()
    for (refmut, resistance) in table.items():
        if refmut < mutations:
            result.append((refmut, resistance))

    return result

def group_resistances(resistances):
    bydrug = defaultdict(list)
    for (mutset, kind) in resistances:
        for degree, drug in zip(kind, ["Oseltamivir", "Zanamivir", "Peramivir", "Laninamivir"]):
            if degree != "NI": # normal inibition
                bydrug[drug].append((mutset, degree))

    return bydrug

def print_report(outpath, bydrug):
    with open(outpath, "w") as file:
        for drug, mutations in bydrug.items():
            print(drug, file=file)
            for mutset, degree in mutations:
                print('\t' + '+'.join(map(str, mutset)), f"({degree})", file=file)
            print("", file=file)

with open(snakemake.input.sub) as file:
    subtype = next(file).strip()

with open(snakemake.input.aln, "rb") as file:
    entries = tools.byte_iterfasta(file)
    reference = next(entries)
    consensus = next(entries)

resistance = get_resistances(reference, consensus, subtype)
bydrug = group_resistances(resistance)

print_report(snakemake.output[0], bydrug)
