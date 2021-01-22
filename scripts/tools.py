import os
import re
import gzip

CODE = {'AAA': 'K', 'CAA': 'Q', 'GAA': 'E', 'TAA': '*',
        'ACA': 'T', 'CCA': 'P', 'GCA': 'A', 'TCA': 'S',
        'AGA': 'R', 'CGA': 'R', 'GGA': 'G', 'TGA': '*',
        'ATA': 'I', 'CTA': 'L', 'GTA': 'V', 'TTA': 'L',
        'AAC': 'N', 'CAC': 'H', 'GAC': 'D', 'TAC': 'Y',
        'ACC': 'T', 'CCC': 'P', 'GCC': 'A', 'TCC': 'S',
        'AGC': 'S', 'CGC': 'R', 'GGC': 'G', 'TGC': 'C',
        'ATC': 'I', 'CTC': 'L', 'GTC': 'V', 'TTC': 'F',
        'AAG': 'K', 'CAG': 'Q', 'GAG': 'E', 'TAG': '*',
        'ACG': 'T', 'CCG': 'P', 'GCG': 'A', 'TCG': 'S',
        'AGG': 'R', 'CGG': 'R', 'GGG': 'G', 'TGG': 'W',
        'ATG': 'M', 'CTG': 'L', 'GTG': 'V', 'TTG': 'L',
        'AAT': 'N', 'CAT': 'H', 'GAT': 'D', 'TAT': 'Y',
        'ACT': 'T', 'CCT': 'P', 'GCT': 'A', 'TCT': 'S',
        'AGT': 'S', 'CGT': 'R', 'GGT': 'G', 'TGT': 'C',
        'ATT': 'I', 'CTT': 'L', 'GTT': 'V', 'TTT': 'F'}

def translate(string, code=CODE):
    "Translates, yielding 'X' if not present."
    if len(string) % 3 != 0:
        raise ValueError("String length must be divisible by 3")

    codons = zip(*[iter(string.upper())]*3)
    return ''.join([code.get("".join(i), 'X') for i in codons])

class Reader:
    """Use this instead of `open` to open files which are either plain text or
    gzipped
    Usage:
    >>> with Reader(file, readmode) as file: # by default textmode
    >>>     print(next(file))
    TEST LINE
    """

    def __init__(self, filename, readmode='r'):
        if readmode not in ('r', 'rt', 'rb'):
            raise ValueError("the Reader cannot write, set mode to 'r' or 'rb'")
        if readmode == 'r':
            self.readmode = 'rt'
        else:
            self.readmode = readmode

        self.filename = filename

        with open(self.filename, 'rb') as f:
            signature = f.peek(8)[:8]

        # Gzipped files begin with the two bytes 0x1F8B
        if tuple(signature[:2]) == (0x1F, 0x8B):
            self.filehandle = gzip.open(self.filename, self.readmode)

        # Else we assume it's a text file.
        else:
            self.filehandle = open(self.filename, self.readmode)

    def close(self):
        self.filehandle.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __iter__(self):
        return self.filehandle

    def __next__(self):
        return next(self.filehandle)

class FastaEntry:
    """One single FASTA entry. Instantiate with string header and bytearray
    sequence."""

    __slots__ = ['header', 'sequence']

    def __init__(self, header, sequence):
        if len(header) > 0 and (header[0] in ('>', '#') or header[0].isspace()):
            raise ValueError('Header cannot begin with #, > or whitespace')
        if '\t' in header:
            raise ValueError('Header cannot contain a tab')

        self.header = header
        self.sequence = sequence.translate(None, b' \t\n\r')

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return '>{}\n{}'.format(self.header, self.sequence.decode())

    def format(self, width=60):
        sixtymers = range(0, len(self.sequence), width)
        spacedseq = '\n'.join([self.sequence[i: i+width].decode() for i in sixtymers])
        return '>{}\n{}'.format(self.header, spacedseq)

    def __getitem__(self, index):
        return self.sequence[index]

    def __repr__(self):
        return '<FastaEntry {}>'.format(self.header)

def byte_iterfasta(filehandle, comment=b'#'):
    """Yields FastaEntries from a binary opened fasta file.
    Usage:
    >>> with Reader('/dir/fasta.fna', 'rb') as filehandle:
    ...     entries = byte_iterfasta(filehandle) # a generator
    Inputs:
        filehandle: Any iterator of binary lines of a FASTA file
        comment: Ignore lines beginning with any whitespace + comment
    Output: Generator of FastaEntry-objects from file
    """

    # Make it work for persistent iterators, e.g. lists
    line_iterator = iter(filehandle)
    # Skip to first header
    try:
        for probeline in line_iterator:
            stripped = probeline.lstrip()
            if stripped.startswith(comment):
                pass

            elif probeline[0:1] == b'>':
                break

            else:
                raise ValueError('First non-comment line is not a Fasta header')

        else: # no break, empty file
            return

    except TypeError:
        errormsg = 'First line does not contain bytes. Are you reading file in binary mode?'
        raise TypeError(errormsg) from None

    header = probeline[1:].rstrip().decode()
    buffer = list()

    # Iterate over lines
    for line in line_iterator:
        if line.startswith(comment):
            pass

        elif line.startswith(b'>'):
            yield FastaEntry(header, bytearray().join(buffer))
            buffer.clear()
            header = line[1:].rstrip().decode()

        else:
            buffer.append(line)

    yield FastaEntry(header, bytearray().join(buffer))

def cat_fasta(outpath, inpaths):
    with open(outpath, "w") as outfile:
        for inpath in inpaths:
            with open(inpath, "rb") as infile:
                for entry in byte_iterfasta(infile):
                    print(entry.format(), file=outfile)

def raise_nonform_fastq(path):
    raise ValueError("FASTQ file {} not correctly formatted".format(filename))

def get_read_pairs(directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError(directory)

    try:
        filenames = sorted(next(os.walk(directory))[2])
    except StopIteration:
        raise FileNotFoundError("Must specify at least one read pair. "
                                "Are you sure you got the directory correct? "
                                f"{os.path.abspath(directory)}")

    pattern = re.compile("(.+)_R([12])(_\d+)?\.fastq\.gz")

    reads = dict()
    for filename in filenames:
        match = pattern.match(filename)
        if match is None:
            raise ValueError("Filename {} does not match pattern \"(.+)_R([12])(_\d+)?\.fastq\.gz\"".format(filename))

        basename, orientation, _ = match.groups()

        if basename not in reads:
            reads[basename] = []
        reads[basename].append(os.path.join(directory, filename))

    singletons = [v for v in reads.values() if len(v) != 2]
    if len(singletons) > 0:
        print("Non-paired read files")
        for group in singletons:
            for file in group:
                print(file)
            print("")

        raise ValueError("Some files not found in pairs")

    return reads

def get_nanopore_reads(directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"{directory} is not a directory.")
    filenames = sorted(next(os.walk(directory))[2])
    if len(filenames) == 0:
        raise ValueError(f"No files found in {directory}")

    pattern = re.compile("([^\.]+)\.fastq\.gz")
    reads = dict()
    for filename in filenames:
        match = pattern.match(filename)
        if match is None:
            raise ValueError(f"Filename {filename} does not match pattern {pattern.pattern}")
        basename = match.groups()[0]
        reads[basename] = os.path.join(directory, filename)
    
    return reads

def fieldsof(filename, lines, n_fields, line_offset):
    for lineno, line in enumerate(map(str.strip, lines)):
        if not line:
            continue

        fields = line.split()
        if len(fields) != n_fields:
            print(fields)
            raise ValueError("Expected {} fields in line {} of {}".format(
                            n_fields, lineno+line_offset, filename))

        yield fields

def get_best_subject(filename):
    with open(filename) as filehandle:
        header = next(filehandle)
        try:
            subject, num, *_ = next(filehandle).strip().split('\t')
        except StopIteration:
            return None, None

    return subject, num

def chain_fasta(paths):
    for path in paths:
        with open(path, "rb") as file:
            for entry in byte_iterfasta(file):
                yield entry

def get_subtype(mappath, subject):
    with open(mappath) as file:
        filtered_lines = filter(str.strip, file)
        splitlines = map(lambda line: line.rstrip().rpartition('\t'), filtered_lines)
        subtype = [t for (s,_,t) in splitlines if s == subject]

    if len(subtype) > 1:
        raise ValueError("Consensus {} has multiple possible subtypes".format(subject))

    if len(subtype) == 0:
        raise ValueError("Consensus {} not found in map".format(subject))

    return subtype[0]

def find_orf(entry, includestop=False):
    """Finds longest ATG xxxx with no STOP pattern in sequence. Includestop to also
    match a terminal STOP codon.
    Returns (pos, orf), with pos being zero-based start pos of orf.
    """

    dna = entry.sequence.decode().upper()
    if dna.translate(str.maketrans("", "", "ATGCWRSYKMBDHVN")):
        raise ValueError("Can only contain ATGCWRSYKMBDHVN")

    codon = "(?:[ATGCWRSYKMBDHVN]{3}(?<!TAG|TAA|TGA))" # codon, not stop
    stop = "(?:TAA|TGA|TAG)"
    start = "ATG"

    if includestop:
        pattern = re.compile(f"(?=({start}{codon}*{stop}))")
    else:
        pattern = re.compile(f"(?=({start}{codon}*))")

    pos, orf = None, ''
    match = max(re.finditer(pattern, dna), key=lambda x: len(x.groups()[0]))
    this_orf = match.groups()[0]
    if len(this_orf) > len(orf):
        orf = this_orf
        pos = match.start()

    orf = FastaEntry(entry.header, bytearray(orf.encode()))
    return pos, orf
