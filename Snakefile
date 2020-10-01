"""
To do:
# IMPROVE REPORT CODE
# Make two distinct output dirs: One containing "refoutput", outputs derived
# only from non-FASTQ files, i.e. the references

Make a gathered resistance_report.txt

* Udtræk Excel ark fra BioNumerics som CSV for en kørsel
* Search beginning of basename in the LIMSnr field
* Get the associated "key"
* For all accepted consensus, for each gene, output a single FASTA file
  with the key as header, and the sequence

  Human er KUN dem med LIMS numre (accession numre) (not those with V)

  * Upload to GISAID
  * Add information to BioNumerics
  *

  fasta header must be key to import to bionumerics
  multiple seqs of same gene can be imported in one go, but not different genes

  file -> import -> fasta from text file -> browse -> pick right "experiment"

  Do not pick FASTA from BioNumerics, instead extract metadata from bionumerics,
  then write a program to extract the data

  Add to bionumerics: SeqCharacterisation (clade), Resistant (yes/no),
  Resistance Mutations (e.g. "K12Y, M19A")


  ------
  Send HA consensus seq for #59 EISN to Ramona
"""

import sys
import os
import itertools as it

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools

# Datadir should look like this:
# ref/
#     H1N1/
#         HA.ref.fna (not including stop)
#         HA.all.fna
#         NA.ref.fna (not including stop)
#         NA.all.fna
#     ...
#     reads/
#         XXXXXX_1.fq.gz
#         XXXXXX_2.fq.gz
#         ...

######################################################
# GLOBAL CONSTANTS
######################################################
READDIR = config["readdir"]
READ_PAIRS = tools.get_read_pairs(os.path.join(SNAKEDIR, READDIR))
BASENAMES = sorted(READ_PAIRS.keys())

REFDIR = os.path.join(SNAKEDIR, "ref")
REFOUTDIR = os.path.join(SNAKEDIR, "refout")
SUBTYPES = set(next(os.walk(REFDIR))[1])
if len(SUBTYPES) < 1:
    raise ValueError("Must specify at least one subtype in data directory")

allsegments = []
for subtype in SUBTYPES:
    filenames = os.listdir(os.path.join(REFDIR, subtype))
    sortedsegments = sorted([fn.partition('.')[0] for fn in filenames if fn.endswith(".ref.fna")])
    allsegments.append((subtype, sortedsegments))

for (subtype, segments) in allsegments[1:]:
    if segments != allsegments[0][1]:
        raise ValueError("Subtype {} does not have the same reference segments as {}".format(subtype, allsegments[0][0]))

SEGMENTS = allsegments[0][1]

CPUS = os.cpu_count()

for subtype in SUBTYPES:
    for segment in SEGMENTS:
        path = REFDIR + "/{}/{}.all.fna".format(subtype, segment)
        if not os.path.exists(path):
            with open(path, "w") as file:
                pass

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

            # Add phylogeny of selected genes
            for phylogene in ["HA"]:
                if phylogene in accepted:
                    inputs.append(f"phylogeny/{basename}/{phylogene}.treefile")

            # Add resistance of selected genes.
            for resistancegene in ["NA"]:
                with open(f"consensus/{basename}/{resistancegene}.subtype") as file:
                    subtype = next(file).strip()

                if resistancegene in accepted and subtype in ("H1N1", "H3N2"):
                    inputs.append(f"consensus/{basename}/{resistancegene}.resistance.txt")

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
# First we do this so we have an "pure subtype" clade.
rule cat_within_subtype:
    input:
        ref=REFDIR + "/{subtype}/{gene}.ref.fna",
        all=REFDIR + "/{subtype}/{gene}.all.fna"
    output: REFOUTDIR + "/{subtype}/{gene}.cat.fna"
    run:
        tools.cat_fasta(output[0], [input.ref, input.all])

rule cat_all_subtypes:
    input: expand(REFOUTDIR + "/{subtype}/{{gene}}.cat.fna", subtype=SUBTYPES),
    output:
        pan=REFOUTDIR + "/pan/{gene,[A-Z0-9]+}.fna",
        map=REFOUTDIR + "/pan/{gene,[A-Z0-9]+}.map.txt"
    run:
        with open(output[0], "w") as panfile, open(output[1], "w") as mapfile:
            for subtype in SUBTYPES:
                path = REFOUTDIR + "/{}/{}.cat.fna".format(subtype, wildcards.gene)
                with open(path, "rb") as infile:
                    for record in tools.byte_iterfasta(infile):
                        print(record.format(), file=panfile)
                        print(record.header, subtype, sep='\t', file=mapfile)

rule index_panfile:
    input: rules.cat_all_subtypes.output.pan
    output:
        comp=REFOUTDIR + "/pan/{gene,[A-Z0-9]+}.comp.b",
        name=REFOUTDIR + "/pan/{gene,[A-Z0-9]+}.name",
        length=REFOUTDIR + "/pan/{gene,[A-Z0-9]+}.length.b",
        seq=REFOUTDIR + "/pan/{gene,[A-Z0-9]+}.seq.b"
    params: REFOUTDIR + "/pan/{gene}"
    log: "log/kma_ref/{gene}.log"
    shell: "kma index -i {input} -o {params} 2> {log}"

rule reference_iqtree:
    input: REFOUTDIR + "/{subtype}/{gene}.cat.aln.fna"
    output: REFOUTDIR + "/{subtype}/{gene}.contree"
    params:
        pre=REFOUTDIR + "/{subtype}/{gene}",
        bootstrap=1000,
        boot_iter=2500,
        model="HKY+G2" # just pick a model to be consistent
    threads: 2
    log: "log/iqtree/{subtype}/{gene}.log"
    shell: "iqtree -s {input} -pre {params.pre} -nt {threads} -m {params.model} "
           "-nm {params.boot_iter} -bb {params.bootstrap} > {log}"

############################
# CONSENSUS PART OF PIPELINE
############################
rule adapterremoval:
    input:
        fw=lambda wildcards: READ_PAIRS[wildcards.basename][0],
        rv=lambda wildcards: READ_PAIRS[wildcards.basename][1],
    output:
        discarded=temp('trim/{basename}/{basename}.discarded.gz'),
        single=temp('trim/{basename}/{basename}.singleton.truncated.gz'),
        fw=temp('trim/{basename}/{basename}.pair1.truncated.gz'),
        rv=temp('trim/{basename}/{basename}.pair2.truncated.gz')
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
        '--minlength 30 --trimns --trimqualities --minquality 20 '
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
        index=rules.index_panfile.output
    output: "aln/{basename}/{gene,[A-Z0-9]+}.spa"
    params:
        db=REFOUTDIR + "/pan/{gene}", # same as index_panfile param
        outbase="aln/{basename}/{gene}"
    threads: 2
    log: "log/aln/{basename}_{gene}.initial.log"
    shell:
        "kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
        "-t {threads} -mmap -Sparse 2> {log}"

def kma_map_index(wc):
    index = tools.get_best_subject(f"aln/{wc.basename}/{wc.gene}.spa")[1]
    return index if index is not None else 1

# Possibly add -ref_fsa to disallow gaps and -dense to disallow insertions
# Align to only the best template found by initial mapping round
rule kma_map:
    input:
        fw='trim/{basename}/{basename}.pair1.truncated.gz',
        rv='trim/{basename}/{basename}.pair2.truncated.gz',
        index=rules.index_panfile.output,
        spa=rules.initial_kma_map.output
    output:
        mat="aln/{basename}/{gene,[A-Z0-9]+}.mat.gz",
        res="aln/{basename}/{gene,[A-Z0-9]+}.res",
        fsa="aln/{basename}/{gene,[A-Z0-9]+}.fsa"
    params:
        db=REFOUTDIR + "/pan/{gene}", # same as index_panfile param
        outbase="aln/{basename}/{gene}",
        refindex=kma_map_index
    log: "log/aln/{basename}_{gene}.log"
    threads: 2
    # This is a run command, because the comamdn cannot be evaluated until
    # the initial_kma_map, so printing it in a shell command will raise an error
    run:
        shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
        "-t {threads} -nf -matrix -mmap -Mt1 {params.refindex} 2> {log}")

rule move_consensus:
    input:
        con="aln/{basename}/{gene}.fsa",
        map=REFOUTDIR + "/pan/{gene}.map.txt"
    output:
        con="consensus/{basename}/{gene}.fna",
        sub="consensus/{basename}/{gene}.subtype"
    run:
        with open(input[0], "rb") as file:
            try:
                consensus = next(tools.byte_iterfasta(file))
                subtype = tools.get_subtype(input.map, consensus.header)
            except StopIteration:
                consensus = tools.FastaEntry("", bytearray())
                subtype = None

        with open(output.sub, "w") as file:
            print(subtype, file=file)

        consensus.header = "{}_{}".format(wildcards.basename, wildcards.gene)
        with open(output.con, "w") as file:
            print(consensus.format(), file=file)

# Here, probably write subtype

rule create_report:
    input:
        con=expand("consensus/{{basename}}/{gene}.fna", gene=SEGMENTS),
        mat=expand("aln/{{basename}}/{gene}.mat.gz", gene=SEGMENTS),
        sub=expand("consensus/{{basename}}/{gene}.subtype", gene=SEGMENTS)
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

rule cat_orfs:
    input:
        consensus="consensus/{basename}/{gene}.orf.faa",
        fnas=expand(REFOUTDIR + "/{subtype}/{{gene}}.cat.orf.faa", subtype=SUBTYPES),
        sub="consensus/{basename}/{gene}.subtype"
    output: "phylogeny/{basename}/{gene,[A-Z0-9]+}.faa"
    run:
        with open(input.sub) as file:
            subtype = next(file).strip()

        tools.cat_fasta(output[0], [f"{REFOUTDIR}/{subtype}/{wildcards.gene}.cat.orf.faa", input.consensus])

def get_iqtree_constraint(wildcards):
    with open(f"consensus/{wildcards.basename}/{wildcards.gene}.subtype") as file:
        subtype = next(file).strip()

    return "{}/{}/{}.contree".format(REFOUTDIR, subtype, wildcards.gene)

rule iqtree:
    input:
        aln="phylogeny/{basename}/{gene}.aln.faa",
        con=get_iqtree_constraint,
        sub="consensus/{basename}/{gene}.subtype"
    output: "phylogeny/{basename}/{gene,[A-Z0-9]+}.treefile"
    params:
        pre="phylogeny/{basename}/{gene}",
        bootstrap=1000,
        boot_iter=2500,
        model="FLU+G2"
    threads: 2
    log: "log/iqtree/phylogeny/{basename}_{gene}.log"
    shell: "iqtree -s {input.aln} -pre {params.pre} -nt {threads} -m {params.model} "
           "-nm {params.boot_iter} -bb {params.bootstrap} -g {input.con} > {log}"

############################
# MUTATIONS PART OF PIPELINE
############################
rule copy_ref:
    input: REFDIR + "/{subtype}/{gene}.ref.fna"
    output: REFOUTDIR + "/{subtype}/{gene}.ref.fna"
    shell: "cp {input} {output}"

rule cat_mutations:
    input:
        con="consensus/{basename}/{gene}.orf.faa",
        refs=expand(REFOUTDIR + "/{subtype}/{{gene}}.ref.orf.faa", subtype=SUBTYPES),
        sub="consensus/{basename}/{gene}.subtype"
    output: "consensus/{basename}/{gene}.cat.faa"
    run:
        with open(input.sub) as file:
            subtype = next(file).strip()
        inpath = f"{REFOUTDIR}/{subtype}/{wildcards.gene}.ref.orf.faa"
        tools.cat_fasta(output[0], [inpath, input.con])

rule mutations:
    input:
        aln="consensus/{basename}/{gene}.cat.aln.faa",
        sub="consensus/{basename}/{gene}.subtype"
    output: "consensus/{basename}/{gene}.resistance.txt"
    script: "scripts/mutations.py"
