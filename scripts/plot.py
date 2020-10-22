import sys
import os
import gzip
import matplotlib.pyplot as plt
import numpy as np

def get_depths(filehandle):
    # Get by-base depth
    depths = list()
    # Try to skip header. If no header exists, no reads were mapped.
    try:
        next(filehandle)
    except StopIteration:
        return [0, 0, 0]

    for fields in filter(None, map(lambda line: line.rstrip().split(), filehandle)):
        # Skip bases that are gap in reference (i.e. insertions)
        if fields[0] == '-':
            continue

        # A,C,G or T
        depths.append(sum(map(int, fields[1:6])))

    return depths

def plot_depths(basepath):
    basename = os.path.basename(basepath)
    genes = set(f.partition('.')[0] for f in os.listdir(basepath) if f.endswith('.mat.gz'))
    figure = plt.figure(dpi=150, figsize=(8, 6))
    for gene in sorted(genes):
        path = os.path.join(basepath, gene + '.mat.gz')
        with gzip.open(path, "rt") as filehandle:
            depths = get_depths(filehandle)

        plt.semilogy(np.linspace(0, 1, len(depths)), depths)

    plt.title(basename, fontsize=14)
    plt.ylim(1, plt.ylim()[1])
    plt.yticks(fontsize=12)
    plt.xlim(-0.01, 1.01)
    plt.legend(sorted(genes))

    return figure

figure = plot_depths(snakemake.input[0])
figure.savefig(snakemake.output[0])
