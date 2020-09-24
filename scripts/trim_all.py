import tools

with open(snakemake.input[0], "rb") as file:
    entries = list(tools.byte_iterfasta(file))

refseq = entries[0].sequence.decode()
leading_gaps = len(refseq) - len(refseq.lstrip('-'))
trailing_gaps = len(refseq) - len(refseq.rstrip('-'))

with open(snakemake.output[0], "w") as file:
    for entry in entries:
        if trailing_gaps == 0:
            newseq = entry.sequence[leading_gaps:]
        else:
            newseq = entry.sequence[leading_gaps:-trailing_gaps]

        entry.sequence = newseq
        print(entry, file=file)
