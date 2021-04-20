import os
import re
import gzip

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

def get_read_pairs(directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Could not locate read dirctory: {os.path.abspath(directory)}")

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

    reads = dict()
    for filename in filenames:
        for ending in [".fq.gz", ".fq", ".fastq.gz", ".fastq"]:
            if filename.endswith(ending):
                reads[filename[:-len(ending)]] = os.path.join(directory, filename)
                break
        else: # no break
            raise ValueError(f"File {filename} is not named like a FASTQ file.")
    
    return reads
