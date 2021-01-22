import sys
import os
import itertools as it

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools


######################################################
# GLOBAL CONSTANTS
######################################################
if "readdir" not in config:
    raise KeyError("You must supply absolute read path: '--config readdir=/path/to/reads'")

if "refset" not in config:
    raise KeyError("You must supply name of reference set: '--config refset=swine'")

READDIR = config["readdir"]
READ_PAIRS = tools.get_read_pairs(READDIR)
BASENAMES = sorted(READ_PAIRS.keys())

REFDIR = os.path.join(SNAKEDIR, "ref", config["refset"])
if not os.path.isdir(REFDIR):
    raise NotADirectoryError(f"Could not find reference directory: '{REFDIR}'")

REFOUTDIR = os.path.join(SNAKEDIR, "refout", config["refset"])

# We have to create these directories either in a rule or outside the DAG to not
# mess up the DAG.
if not os.path.isdir(os.path.join(SNAKEDIR, "refout")):
    os.mkdir(os.path.join(SNAKEDIR, "refout"))

if not os.path.isdir(os.path.join(SNAKEDIR, "refout", config["refset"])):
    os.mkdir(os.path.join(SNAKEDIR, "refout", config["refset"]))

SEGMENTS = sorted([fn.rpartition('.')[0] for fn in os.listdir(REFDIR)
    if fn.endswith(".fna")])

######################################################
# Start of pipeline
######################################################
ruleorder: translate > trim_alignment > mafft

# We add the checkpoiint output here just to make sure the checkpoint is included
# in the DAG
def done_input(wildcards):
    # Add report
    inputs = ["report.txt"]
    checkpoints.cat_reports.get()

    for basename in BASENAMES:
        # Add FASTQC report
        for pair in (1, 2):
            inputs.append(f"fastqc/{basename}/{basename}.pair{pair}.truncated_fastqc.html")
        
        # Add individual report
        inputs.append(f"consensus/{basename}/report.txt")

        # Add depth plot
        inputs.append(f"depths/{basename}.pdf")
        with open(f"consensus/{basename}/accepted.txt") as accepted_file:
            accepted = set(filter(None, map(str.strip, accepted_file)))

            # Add phylogeny of selected segments
            """
            for phylosegment in ["HA"]:
                if phylosegment in accepted:
                    inputs.append(f"phylogeny/{basename}/{phylosegment}.treefile")

            # Add resistance of selected segments.
            for resistancesegment in ["NA"]:
                with open(f"consensus/{basename}/{resistancesegment}.subtype") as file:
                    subtype = next(file).strip()

                if resistancesegment in accepted and subtype in ("H1N1", "H3N2"):
                    inputs.append(f"consensus/{basename}/{resistancesegment}.resistance.txt")
            """

    return inputs

rule all:
    input: done_input
    output: touch("done")

############################
# Alignment module
############################
rule mafft:
    input: "{dir}/{base}.{ext}"
    output: "{dir}/{base}.aln.{ext,faa|fna}"
    log: "log/mafft/{dir}/{base}.{ext}.log"
    shell: "mafft --auto {input} > {output} 2> {log}"

rule trim_alignment:
    input: "{dir}/{base}.aln.{ext}"
    output: "{dir}/{base}.aln.trimmed.{ext,faa|fna}"
    script: "scripts/trim_all.py"

rule translate:
    input: "{dir}/{base}.fna"
    output: "{dir}/{base}.faa"
    run:
        with tools.Reader(input[0], "rb") as file:
            entries = list(tools.byte_iterfasta(file))

        with open(output[0], "w") as file:
            for entry in entries:
                translated = tools.translate(entry.sequence.decode())
                print(">{}\n{}".format(entry.header, translated), file=file)

rule extract_orf:
    input: "{dir}/{base}.fna"
    output: "{dir}/{base}.orf.fna"
    run:
        with open(input[0], "rb") as file, open(output[0], "w") as outfile:
            for entry in tools.byte_iterfasta(file):
                pos, orf = tools.find_orf(entry)
                print(orf.format(), file=outfile)

#################################
# REFERENCE-ONLY PART OF PIPELINE
#################################
# We use k=14 because we sometimes have very small reads.
rule index_ref:
    input: REFDIR + "/{segment}.fna"
    output:
        comp=REFOUTDIR + "/{segment,[a-zA-Z0-9]+}.comp.b",
        name=REFOUTDIR + "/{segment,[a-zA-Z0-9]+}.name",
        length=REFOUTDIR + "/{segment,[a-zA-Z0-9]+}.length.b",
        seq=REFOUTDIR + "/{segment,[a-zA-Z0-9]+}.seq.b"
    params:
        outpath=REFOUTDIR + "/{segment}",
    log: "log/kma_ref/{segment}.log"
    shell: "kma index -k 14 -i {input} -o {params.outpath} 2> {log}"

""" TODO: Create ref tree for each "close" group of refs.
rule reference_iqtree:
    input: REFOUTDIR + "/{subtype}/{segment}.cat.aln.fna"
    output: REFOUTDIR + "/{subtype}/{segment}.contree"
    params:
        pre=REFOUTDIR + "/{subtype}/{segment}",
        bootstrap=1000,
        boot_iter=2500,
        model="HKY+G2" # just pick a model to be consistent
    threads: 2
    log: "log/iqtree/{subtype}/{segment}.log"
    shell: "iqtree -s {input} -pre {params.pre} -nt {threads} -m {params.model} "
           "-nm {params.boot_iter} -bb {params.bootstrap} > {log}"
"""

############################
# CONSENSUS PART OF PIPELINE
############################
rule adapterremoval:
    input:
        fw=lambda wildcards: READ_PAIRS[wildcards.basename][0],
        rv=lambda wildcards: READ_PAIRS[wildcards.basename][1],
    output:
        discarded='trim/{basename}/{basename}.discarded.gz',
        single='trim/{basename}/{basename}.singleton.truncated.gz',
        fw='trim/{basename}/{basename}.pair1.truncated.gz',
        rv='trim/{basename}/{basename}.pair2.truncated.gz'
    log: 'log/trim/{basename}.log'
    params:
        basename="trim/{basename}/{basename}"
    threads: 2
    shell:
        'AdapterRemoval '
        # Input files
        '--file1 {input.fw} --file2 {input.rv} '
        # Output files
        '--basename {params.basename} '
        "2> {log} "
        # Other parameters:
        '--minlength 20 --trimns --trimqualities --minquality 20 '
        '--qualitybase 33 --gzip --trimwindows 5 --threads {threads}'

rule fastqc:
    input: [rules.adapterremoval.output.fw, rules.adapterremoval.output.rv]
    output:
        fw='fastqc/{basename}/{basename}.pair1.truncated_fastqc.html',
        rv='fastqc/{basename}/{basename}.pair2.truncated_fastqc.html',
    log: 'log/fastqc/{basename}.log'
    threads: 2
    # Annoyingly, it prints "analysis complete" to stdout
    shell: 'fastqc -t {threads} {input} -o fastqc/{wildcards.basename} 2> {log} > /dev/null'

# Do this to get the best template, -Sparse option is designed for this.
rule initial_kma_map:
    input:
        fw='trim/{basename}/{basename}.pair1.truncated.gz',
        rv='trim/{basename}/{basename}.pair2.truncated.gz',
        # We artificially add FASTQC here to move FASTQC up the dependency graph,
        # such that it is completed before the checkpoint. Else, Snakemake will remove
        # the trimmed files before FASTQC, since it can't see across the checkpoint
        # and will believe there is no more use of the trimmed files.
        fastqc = rules.fastqc.output,
        index=rules.index_ref.output
    output: "aln/{basename}/{segment,[A-Z0-9]+}.spa"
    params:
        db=REFOUTDIR + "/{segment}", # same as index_reffile param
        outbase="aln/{basename}/{segment}"
    threads: 2
    log: "log/aln/{basename}_{segment}.initial.log"
    shell:
        "kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
        "-t {threads} -Sparse 2> {log}"

def kma_map_index(wc):
    index = tools.get_best_subject(f"aln/{wc.basename}/{wc.segment}.spa")[1]
    return index if index is not None else 1

# Possibly add -ref_fsa to disallow gaps and -dense to disallow insertions
# Align to only the best template found by initial mapping round
rule kma_map:
    input:
        fw='trim/{basename}/{basename}.pair1.truncated.gz',
        rv='trim/{basename}/{basename}.pair2.truncated.gz',
        index=rules.index_ref.output,
        spa=rules.initial_kma_map.output
    output:
        mat="aln/{basename}/{segment,[A-Z0-9]+}.mat.gz",
        res="aln/{basename}/{segment,[A-Z0-9]+}.res",
        fsa="aln/{basename}/{segment,[A-Z0-9]+}.fsa"
    params:
        db=REFOUTDIR + "/{segment}", # same as index_reffile param
        outbase="aln/{basename}/{segment}",
        refindex=kma_map_index
    log: "log/aln/{basename}_{segment}.log"
    threads: 2
    # This is a run command, because the comamdn cannot be evaluated until
    # the initial_kma_map, so printing it in a shell command will raise an error
    # If gapopen is default, it will lead to totally absurd alignments. Perhaps
    # setting open/ext/m/mm be -5/-1/1/-3 will lead to better results still?
    run:
        shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
        "-t {threads} -k 16 -gapopen -5 -nf -matrix -Mt1 {params.refindex} 2> {log}")

rule move_consensus:
    input: "aln/{basename}/{segment}.fsa"
    output: "consensus/{basename}/{segment}.untrimmed.fna"
    run:
        with open(input[0], "rb") as file:
            try:
                consensus = next(tools.byte_iterfasta(file))
            except StopIteration:
                consensus = tools.FastaEntry("", bytearray())

        consensus.header = "{}_{}".format(wildcards.basename, wildcards.segment)
        with open(output[0], "w") as file:
            print(consensus.format(), file=file)

# In our current lab setup, we use primers to amplify our influeza segments. But these do not have the
# proper sequence
rule remove_primers:
    input:
        con=rules.move_consensus.output,
        primers=f"{SNAKEDIR}/ref/primers.fna"
    output: "consensus/{basename}/{segment}.fna"
    log: "log/consensus/remove_primers_{basename}_{segment}.txt"
    params:
        scriptpath=f"{SNAKEDIR}/scripts/trim_consensus.jl",
        minmatches=6
    run:
        shell("julia {params.scriptpath} {input.primers} {input.con} {output} {params.minmatches} > {log}")

rule create_report:
    input:
        con=expand("consensus/{{basename}}/{segment}.fna", segment=SEGMENTS),
        mat=expand("aln/{{basename}}/{segment}.mat.gz", segment=SEGMENTS)
    output:
        rep="consensus/{basename}/report.txt",
        acc="consensus/{basename}/accepted.txt"
    script: "scripts/report.py"

checkpoint cat_reports:
    input: expand("consensus/{basename}/report.txt", basename=BASENAMES)
    output: "report.txt"
    run:
        with open(output[0], "w") as outfile:
            for basename in BASENAMES:
                with open(f"consensus/{basename}/report.txt") as infile:
                    print(basename, file=outfile)
                    for line in infile:
                        print('\t', line, sep='', end='', file=outfile)

rule plot_depths:
    input: "aln/{basename}"
    output: "depths/{basename}.pdf"
    script: "scripts/plot.py"

############################
# IQTREE PART OF PIPELINE
############################
# We probably only want to run IQ tree on the closest subsequences TBH.
# Later maybe I can just pick out the fitting subtype and run IQTREE on those
# Actually maybe even just IQTREE on the reference then use that as a guide tree

"""
rule cat_orfs:
    input:
        consensus="consensus/{basename}/{segment}.orf.faa",
        fnas=expand(REFOUTDIR + "/{subtype}/{{segment}}.cat.orf.faa", subtype=SUBTYPES),
        sub="consensus/{basename}/{segment}.subtype"
    output: "phylogeny/{basename}/{segment,[A-Z0-9]+}.faa"
    run:
        with open(input.sub) as file:
            subtype = next(file).strip()

        tools.cat_fasta(output[0], [f"{REFOUTDIR}/{subtype}/{wildcards.segment}.cat.orf.faa", input.consensus])

def get_iqtree_constraint(wildcards):
    with open(f"consensus/{wildcards.basename}/{wildcards.segment}.subtype") as file:
        subtype = next(file).strip()

    return "{}/{}/{}.contree".format(REFOUTDIR, subtype, wildcards.segment)

rule iqtree:
    input:
        aln="phylogeny/{basename}/{segment}.aln.faa",
        con=get_iqtree_constraint,
        sub="consensus/{basename}/{segment}.subtype"
    output: "phylogeny/{basename}/{segment,[A-Z0-9]+}.treefile"
    params:
        pre="phylogeny/{basename}/{segment}",
        bootstrap=1000,
        boot_iter=2500,
        model="FLU+G2"
    threads: 2
    log: "log/iqtree/phylogeny/{basename}_{segment}.log"
    shell: "iqtree -s {input.aln} -pre {params.pre} -nt {threads} -m {params.model} "
           "-nm {params.boot_iter} -bb {params.bootstrap} -g {input.con} > {log}"

############################
# MUTATIONS PART OF PIPELINE
############################
rule copy_ref:
    input: REFDIR + "/{subtype}/{segment}.ref.fna"
    output: REFOUTDIR + "/{subtype}/{segment}.ref.fna"
    shell: "cp {input} {output}"

rule cat_mutations:
    input:
        con="consensus/{basename}/{segment}.orf.faa",
        refs=expand(REFOUTDIR + "/{subtype}/{{segment}}.ref.orf.faa", subtype=SUBTYPES),
        sub="consensus/{basename}/{segment}.subtype"
    output: "consensus/{basename}/{segment}.cat.faa"
    run:
        with open(input.sub) as file:
            subtype = next(file).strip()
        inpath = f"{REFOUTDIR}/{subtype}/{wildcards.segment}.ref.orf.faa"
        tools.cat_fasta(output[0], [inpath, input.con])

rule mutations:
    input:
        aln="consensus/{basename}/{segment}.cat.aln.faa",
        sub="consensus/{basename}/{segment}.subtype"
    output: "consensus/{basename}/{segment}.resistance.txt"
    script: "scripts/mutations.py"
"""