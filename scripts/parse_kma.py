import tools
import itertools as it

def parse_matrix(filename, filehandle, subject):
    # First find the correct matrix
    for line in filehandle:
        if line[1:-1] == subject:
            break
    else:
        raise ValueError("Subject {} not found in matrix {}".format(subject, filename))

    lines = it.takewhile(lambda line: line.strip() and line[0] != '#', filehandle)
    matrix = list()
    for fields in tools.fieldsof(filename, lines, 7, 1):
        matrix.append([fields[0]] + [int(i) for i in fields[1:]])

    return matrix

def reconstruct_sequence(matrix):
    sequence = list()
    for position in matrix:
        base = None
        refbase = position[0]

        # If depth < 5, just choose reference base
        counts = position[1:6]
        maxval = max(counts)
        if maxval < 5:
            base = refbase
        elif sorted(counts)[-2] > maxval * 0.5:
            base = refbase
        else:
            # Choose the base with most reads. This is a homemade argmax function
            for n, bs in zip(counts, "ACGTN"):
                if n == maxval:
                    base = bs
                    break

            # Should never happen
            else:
                assert False

        sequence.append(base)

    return ''.join(sequence)

def translate_and_verify(subject, sequence):
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
    subjects = tools.get_possible_subjects(snakemake.input[1], file)

    # We should only have one subject here
    if len(subjects) > 1:
        raise ValueError("More than 1 subject in {}".format(snakemake.input[1]))

    subject = subjects[0]

with tools.Reader(snakemake.input[0]) as file:
    matrix = parse_matrix(snakemake.input[0], file, subject)
    sequence = reconstruct_sequence(matrix)

with open(snakemake.output[0], "w") as file:
    print(">{}\n{}".format(snakemake.wildcards.basename, sequence), file=file)

with open(snakemake.output[1], "w") as file:
    print(">{}\n{}".format(snakemake.wildcards.basename, translate_and_verify(subject, sequence)), file=file)
