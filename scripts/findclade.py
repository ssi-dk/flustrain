import tools
import re

class Mutation:
    """
    Pos is stored zero-indexed such that the start methionine is pos=0
    """
    regex = re.compile("^([A-Z]?)(-?\d+)([A-Z]?)$")
    __slots__ = ['ref', 'pos', 'alt']

    def __init__(self, ref, pos, alt):
        self.ref = ref
        self.pos = pos
        self.alt = alt

        if ref is None and alt is None:
            raise ValueError("Both ref and alt cannot be None")

    @classmethod
    def fromstring(cls, st, offset=0):
        match = cls.regex.match(st)
        if match is None:
            raise ValueError("String '{}' is not formatted as mutation".format(st))
        ref, pos, alt = match.groups()
        ref = None if ref == "" else ref
        alt = None if alt == "" else alt
        pos = int(pos) - 1 + offset
        return cls(ref, pos, alt)

    def __repr__(self):
        a = "" if self.ref is None else self.ref
        b = "" if self.alt is None else self.alt
        return "{}{}{}".format(a, self.pos+1, b)

    def __eq__(self, other):
        return self.ref == other.ref and self.pos == other.pos and self.alt == other.alt

    def __hash__(self):
        return hash(self.ref) ^ hash(self.pos) ^ hash(self.alt)

class Clade:
    __slots__ = ['name', 'definition', 'subclades']

    def __init__(self, name, definition, parent=None):
        definition = set(definition)

        self.name = name

        if parent is None:
            self.definition = definition
        else:
            self.definition = definition | parent.definition

        self.subclades = []

        if parent is not None:
            parent.add(self)

    def intersect(self, mutationset):
        return self.definition & mutationset

    def matches(self, mutationset):
        return self.definition.issubset(set(mutationset))

    def add(self, clade):
        if not isinstance(clade, Clade):
            raise ValueError("Can only add Clades to a Clade")
        self.subclades.append(clade)

# HA2 offset in H1N1 is 344, e.g. HA2 X1Y is HA X345Y
H1N1 = {"H1": ([], {
    "6B.1": (["S162N", "I216T"], {
        "6B.1A": (["S74R", "S164T", "I295V"], {
            "6B.1A1": (["S183P"], {}),
            "6B.1A2": (["S183P", "L233I", "HA2 V193A"], {}),
            "6B.1A3": (["T120A", "S183P"], {}),
            "6B.1A4": (["N129D", "A144E", "S183P"], {}),
            "6B.1A5": (["S183P", "N260D"], {
                "6B.1A5A": (["N129D", "T185A"], {}),
                "6B.1A5B": (["E235D", "HA2 V193A"], {}),
            }),
            "6B.1A6": (["T120A", "S183P"], {}),
            "6B.1A7": (["K302T", "HA2 I77M", "HA2 N169S", "HA2 E179D"], {})
        })
    })
})
}

# Note: WHO made a mistake and accidentally listed N145S as N144S
H3N2 = {"H3": ([], {
    "3C.2a": (["L3I", "N145S", "F159Y", "K160T", "N225D", "Q311H", "HA2 D160N"], {
        "3C.2a1": (["N171K", "HA2 I77V", "HA2 G155E"], {
            "3C.2a1a": (["T135K", "HA2 G150E"], {}),
            "3C.2a1b": (["E62G", "R142G", "H311Q"], {}),
        }),
        "3C.2a2": (["T131K", "R142K", "R261Q",], {}),
        "3C.2a3": (["N121K", "S144K"], {}),
        "3C.2a4": (["N31S", "D53N", "S144R", "N171K", "I192T", "Q197H"], {})
    }),
    "3C.3a": (["A138S", "F159S", "N225D", "K326R"], {}),
})
}

# Parse the mutatitions
# This function takes (name: ([], {})) tuple
def parse_clade(parent, pair, offset):
    name, (definition, subclades) = pair
    mutations = []
    for st in definition:
        if st.startswith("HA2"):
            thisoffset = 344
            st = st.partition(" ")[2]
        else:
            thisoffset = offset

        mutation = Mutation.fromstring(st, thisoffset)
        mutations.append(mutation)

    clade = Clade(name, mutations, parent)
    for pair in subclades.items():
        parse_clade(clade, pair, offset)

    return clade

H1_clade = parse_clade(None, next(iter(H1N1.items())), 0)
H3_clade = parse_clade(None, next(iter(H3N2.items())), 0)

# Logic to find the right clade
def findclade(clade, mutations):
    mutset = set(mutations)
    result = clade
    subclades = [c for c in clade.subclades if c.matches(mutations)]

    # Strangely, one subclade A may have a set of mutations which is a subset
    # of clade B, while A not being annotated as a subclade of B. My program chooses
    # the subclade with most mutations in that case.
    if len(subclades) > 1:
        subclades.sort(key=lambda x: len(x.definition), reverse=True)
        if all(sc.definition.issubset(subclades[0].definition) for sc in subclades[1:]):
            subclades = [subclades[0]]

    if len(subclades) == 1:
        return findclade(subclades[0], mutations)
    else:
        return result

if __name__ == '__main__':
    with open(snakemake.input.subtype) as file:
        subtype = next(file).strip()

    offset = tools.get_offset(subtype)
    print(snakemake.input.subtype, subtype, offset)

    mutations = set()
    with open(snakemake.input.mutations) as file:
        for line in file:
            mutations.add(Mutation.fromstring(line.strip(), offset))

    if subtype.startswith("H1"):
        clade = findclade(H1_clade, mutations).name
    elif subtype.startswith("H3"):
        clade = findclade(H3_clade, mutations).name
    else:
        clade = "unknown clade"

    with open(snakemake.output[0], "w") as file:
        print(clade, file=file)
