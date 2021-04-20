# Directory structure of ref/trees is:

# . refset [e.g. human / swine / avian]
# ├── segment1 [segment name: e.g. HA]
# │   ├── subtype1.fna [arbitrary: e.g. H1sw.fna]
# │   └── subtype2.fna [arbitrary: e.g. H2.fna] 
# └── segment2 [e.g. NA]
#     ├── subtype1.fna [arbitrary: e.g. N1.fna]
#     └── subtype2.fna [arbitrary: e.g. N2sw.fna]

# Each segment is processed independently.
# * First, the closest subtype is found by roughly comparing against the sequences
# in each subtype.
# * Then, the segment is added to the subtype fasta to obtain the "cat"
# * Then, the cat is aligned and trimmed
# * Then, the trimmed is run through phylogeny

import sys
import os

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools

######################################################
# GLOBAL CONSTANTS
######################################################
JULIA_COMMAND = f"julia --startup-file=no --project={SNAKEDIR}"

if "refset" not in config or config["refset"] not in ["human", "swine", "avian"]:
    raise KeyError("You must supply reference set: '--config refset=[\"human/swine/avian\"]'")

REFOUTDIR = os.path.join(SNAKEDIR, "refout")

# We have to create this directories either in a rule or outside the DAG to not
# mess up the DAG.
if not os.path.isdir(REFOUTDIR):
    os.mkdir(REFOUTDIR)

REFTREEDIR = os.path.join(REFOUTDIR, config["refset"])
if not os.path.isdir(REFTREEDIR):
    os.mkdir(REFTREEDIR)

REFSEQDIR = os.path.join(SNAKEDIR, "ref", "trees", config["refset"])
if not os.path.isdir(REFSEQDIR):
    raise NotADirectoryError(f"Reference tree sequence directory {REFSEQDIR} not found")

TREESEGMENTS = sorted(os.listdir(REFSEQDIR))

SUBTYPES = dict()
for segment in TREESEGMENTS:
    SUBTYPES[segment] = [i[:-4] for i in os.listdir(os.path.join(REFSEQDIR, segment)) if i.endswith(".fna")]

if "consensus" not in config:
    raise KeyError("You must supply path to consensus directory: '--config consensus=/path/to/consensus'")

if not os.path.isdir(config["consensus"]):
    raise NotADirectoryError(f'Consensus directory {config["consensus"]} not found')

BASENAMES = sorted(os.listdir(config["consensus"]))

PASSED = dict()
for basename in BASENAMES:
    PASSED[basename] = []

    filename = os.path.join(config["consensus"], basename, "curated.fna")
    if not os.path.isfile(filename):
        continue

    with open(filename, "rb") as file:
        for record in tools.byte_iterfasta(file):
            segment = record.header.split('_')[-1]
            if segment in TREESEGMENTS:
                PASSED[basename].append(segment)

#################################
# ALL RULE
#################################

def done_input(wildcards):
    inputs = []
    for (basename, segments) in PASSED.items():
        for segment in segments:
            inputs.append(f"phylogeny/{basename}/{segment}.iqtree.treefile")
    
    return inputs

rule all:
    input: done_input

#################################
# REFERENCE-ONLY PART OF PIPELINE
#################################

# This can be changed from "linsi" to "mafft" if it's too slow
rule guide_tree_alignment:
    input: REFSEQDIR + "/{segment}/{subtype}.fna"
    output: REFTREEDIR + "/{segment}/{subtype}.aln.fna"
    threads: 2
    log: "log/ref_phylogeny/{segment}_{subtype}.aln.log"
    shell: "linsi --quiet --thread {threads} {input} > {output}"

rule guide_tree_trimal:
    input: rules.guide_tree_alignment.output
    output: REFTREEDIR + "/{segment}/{subtype}.aln.trim.fna"
    shell: "trimal -gt 0.9 -cons 60 -keepheader -keepseqs -in {input} -fasta -out {output}"

rule guide_tree_iqtree:
    input: rules.guide_tree_trimal.output
    output: REFTREEDIR + "/{segment}/{subtype}.treefile"
    log: "log/ref_phylogeny/{segment}_{subtype}.iqtree.log"
    threads: 2
    params:
        pre=REFTREEDIR + "/{segment}/{subtype}",
        model="HKY+G2",
    shell: "iqtree -T {threads} --quiet -s {input} -pre {params.pre} -m {params.model}"

#################################
# PHYLOGENY PART OF PIPELINE
#################################

# This rule extracts the relevant sequences from all curated consensus files,
# then finds the best subtype for each.
rule get_subtype:
    input: expand("consensus/{basename}/curated.fna", basename=BASENAMES)
    output:
        subtype=expand("phylogeny/{basename}/{segment}_subtype.txt", basename=BASENAMES, segment=PASSED[basename]),
        consensus=expand("phylogeny/{basename}/{segment}.fna", basename=BASENAMES, segment=PASSED[basename])
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/get_subtype.jl",
        consensus=config["consensus"],
        refseqdir=REFSEQDIR
    shell: "{params.juliacmd} {params.scriptpath} {params.refseqdir} {params.consensus} phylogeny"

def get_subtype(basename, segment):
    with open(f"phylogeny/{basename}/{segment}_subtype.txt") as file:
        lines = list(filter(None, map(str.strip, file)))
        if len(lines) != 1:
            raise ValueError(f"phylogeny/{basename}/{segment}_subtype.txt doesnt have one subtype")
        return lines[0]

rule align_consensus_to_ref:
    input:
        consensus=rules.get_subtype.output.consensus,
        refaln=expand(REFTREEDIR + "/{{segment}}/{subtype}.aln.trim.fna", subtype=SUBTYPES[segment])
    output: temp("phylogeny/{basename}/{segment}.cat.aln.fna")
    log: "log/mafft/phylogeny/{basename}_{segment}.log"
    run:
        subtype = get_subtype(wildcards.basename, wildcards.segment)
        refalnpath = REFTREEDIR + f"/{wildcards.segment}/{subtype}.aln.trim.fna"
        shell("mafft --add phylogeny/{wildcards.basename}/{segment}.fna --keeplength {refalnpath} > {output} 2> {log}")

def all_guidetrees(wildcards):
    return [REFTREEDIR + f"/{wildcards.segment}/{subtype}.treefile" for subtype in SUBTYPES[wildcards.segment]]

def best_guidetree(wildcards):
    subtype = get_subtype(wildcards.basename, wildcards.segment)
    return REFTREEDIR + f"/{wildcards.segment}/{subtype}.treefile"

rule iqtree:
    input:
        consensus=rules.align_consensus_to_ref.output,
        # We need to create all guide trees because the correct one is decided at
        # runtime.
        guidetrees=all_guidetrees
    output: "phylogeny/{basename}/{segment}.iqtree.treefile"
    params:
        guidetree=best_guidetree,
        pre="phylogeny/{basename}/{segment}.iqtree",
        model="HKY+G2"
    threads: 2
    log: "log/iqtree/phylogeny/{basename}_{segment}.log"
    run:
        shell("iqtree -s {input[0]} -pre {params.pre} -nt {threads} -m {params.model} "
            "-g {params.guidetree} > {log}")