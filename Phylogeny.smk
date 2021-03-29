import sys
import os

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools

######################################################
# GLOBAL CONSTANTS
######################################################
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

TREESEGMENTS = sorted([fn.rpartition('.')[0] for fn in os.listdir(REFSEQDIR) if fn.endswith(".fna")])

if "consensus" not in config:
    raise KeyError("You must supply path to consensus directory: '--config consensus=/path/to/consensus'")

if not os.path.isdir(config["consensus"]):
    raise NotADirectoryError(f'Consensus directory {config["consensus"]} not found')

BASENAMES = sorted(os.listdir(config["consensus"]))

#################################
# ALL RULE
#################################

def done_input(wildcards):
    inputs = []
    for basename in BASENAMES:
        filename = os.path.join(config["consensus"], basename, "curated.fna")
        if not os.path.isfile(filename):
            continue

        with open(filename, "rb") as file:
            for record in tools.byte_iterfasta(file):
                segment = record.header.split('_')[-1]
                if segment in TREESEGMENTS:
                    inputs.append(f"phylogeny/{basename}/{segment}.iqtree.treefile")
    
    return inputs

rule all:
    input: done_input

#################################
# REFERENCE-ONLY PART OF PIPELINE
#################################

# This can be changed from "linsi" to "mafft" if it's too slow
rule guide_tree_alignment:
    input: REFSEQDIR + "/{segment}.fna"
    output: REFTREEDIR + "/{segment}.aln.fna"
    threads: 2
    log: "log/ref_phylogeny/{segment}.aln.log"
    shell: "linsi --quiet --thread {threads} {input} > {output}"

rule guide_tree_trimal:
    input: rules.guide_tree_alignment.output
    output: REFTREEDIR + "/{segment}.aln.trim.fna"
    shell: "trimal -gt 0.9 -cons 60 -keepheader -keepseqs -in {input} -fasta -out {output}"

rule guide_tree_iqtree:
    input: rules.guide_tree_trimal.output
    output: REFTREEDIR + "/{segment}.treefile"
    log: "log/ref_phylogeny/{segment}.iqtree.log"
    threads: 2
    params:
        pre=REFTREEDIR + "/{segment}",
        model="HKY+G2",
    shell: "iqtree -T {threads} --quiet -s {input} -pre {params.pre} -m {params.model}"

#################################
# PHYLOGENY PART OF PIPELINE
#################################

rule extract_phylogeny_seq:
    input: "consensus/{basename}/curated.fna"
    output: temp("phylogeny/{basename}/{segment}.fna")
    run:
        theentry = None
        with tools.Reader(input[0], "rb") as file:
            for entry in tools.byte_iterfasta(file):
                if entry.header == f"{wildcards.basename}_{wildcards.segment}":
                    theentry = entry
                    break
        
        if theentry is None:
            raise ValueError(f"For {basename}, tree seq not found in curated.fna: {wildcards.segment}")

        with open(output[0], "w") as file:
            print(entry.format(), file=file)

rule align_consensus_to_ref:
    input:
        consensus="phylogeny/{basename}/{segment}.fna",
        refaln=REFTREEDIR + "/{segment}.aln.trim.fna"
    output: temp("phylogeny/{basename}/{segment}.cat.aln.trim.fna")
    log: "log/mafft/phylogeny/{basename}_{segment}.log"
    shell: "mafft --add {input.consensus} {input.refaln} > {output} 2> {log}"

rule iqtree:
    input:
        consensus=rules.align_consensus_to_ref.output,
        guidetree=rules.guide_tree_iqtree.output
    output: "phylogeny/{basename}/{segment,[A-Z0-9]+}.iqtree.treefile"
    params:
        pre="phylogeny/{basename}/{segment}.iqtree",
        model="HKY+G2"
    threads: 2
    log: "log/iqtree/phylogeny/{basename}_{segment}.log"
    shell: "iqtree -s {input[0]} -pre {params.pre} -nt {threads} -m {params.model} "
           "-g {input[1]} > {log}"