import sys
import os
import itertools as it

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools

# We only need this because Julia 1.6 segfaults with PkgCompiler at the moment
JULIA_PATH = "/Applications/Julia-1.5.app/Contents/Resources/julia/bin/julia"
SYSIMG_PATH = os.path.join(SNAKEDIR, "scripts", "sysimg", "sysimg.so")
JULIA_COMMAND = f"{JULIA_PATH} --startup-file=no --project={SNAKEDIR} -J {SYSIMG_PATH}"

######################################################
# GLOBAL CONSTANTS
######################################################
if "readdir" not in config:
    raise KeyError("You must supply absolute read path: '--config readdir=/path/to/reads'")

if "refset" not in config:
    raise KeyError("You must supply name of reference set: '--config refset=swine'")

if config["refset"] not in ["human", "swine", "avian"]:
    raise KeyError("refset must be 'human', 'swine' or 'avian'")

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

# Only relevant for the human reference set with its fixed subtypes
SUBTYPES = ["H1N1", "H3N2", "Victoria", "Yamagata"]

######################################################
# Start of pipeline
######################################################
ruleorder: translate > trim_alignment > mafft

# We add the checkpoiint output here just to make sure the checkpoint is included
# in the DAG
def done_input(wildcards):
    # Add report and the commit
    inputs = ["report.txt"]
    checkpoints.cat_reports.get()

    for basename in BASENAMES:
        # Add trimmed FASTQ
        for direction in ["fw", "rv"]:
            inputs.append(f"trim/{basename}/{direction}.fq.gz")

        # Add individual report
        inputs.append(f"consensus/{basename}/report.txt")

        # Add depth plot
        inputs.append(f"depths/{basename}.pdf")
        with open(f"consensus/{basename}/accepted.txt") as accepted_file:
            accepted = set(filter(None, map(str.strip, accepted_file)))

        # Add phylogeny - only for human viruses for which we have
        # well defined subtypes. 
        if config["refset"] == "human":
            inputs.append(f"phylogeny/{basename}/HA.treefile")

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
    output: "commit.txt"
    params: SNAKEDIR
    shell: "git -C {params} rev-parse --short HEAD > {output}"

############################
# Alignment module
############################
rule mafft:
    input: "{dir}/{base}.{ext}"
    output: "{dir}/{base}.aln.{ext,faa|fna}"
    threads: 4
    shell: "mafft --thread {threads} --auto {input} > {output}"

rule trim_alignment:
    input: "{dir}/{base}.aln.{ext}"
    output: "{dir}/{base}.aln.trimmed.{ext,faa|fna}"
    shell: "trimal -gt 0.9 -cons 60 -keepheader -keepseqs -in {input} -fasta -out {output}"

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
# Since we don't re-calculate this index for every run, we can't toggle
# the k value depending on the read length config.
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
    
    # The pipeline is very sensitive to the value of k here.
    # Too low means the mapping is excruciatingly slow,
    # too high results in very poor mapping quality.
    shell: "kma index -k 12 -i {input} -o {params.outpath} 2> {log}"


rule split_ref_subtype:
    input: REFDIR + "/{segment}.fna"
    output: expand(REFOUTDIR + "/{{segment}}_{subtype}.fna", subtype=SUBTYPES)
    run:
        with open(input[0], "rb") as infile:
            records = list(tools.byte_iterfasta(infile))
        
        for subtype in SUBTYPES:
            recs = [r for r in records if r.header.split('|')[1] == subtype]
            with open(REFOUTDIR + f"/{wildcards.segment}_{subtype}.fna", "w") as file:
                for record in recs:
                    print(record.format(), file=file)

rule reference_iqtree:
    input: REFOUTDIR + "/{segment}_{subtype}.aln.trimmed.fna"
    output: REFOUTDIR + "/{segment}_{subtype}.treefile"
    params:
        pre=REFOUTDIR + "/{segment}_{subtype}",
        bootstrap=1000,
        boot_iter=2500,
        model="HKY+G2" # just pick a model to be consistent
    threads: 2
    log: "log/iqtree/{subtype}/{segment}.log"
    shell: "iqtree -s {input} -pre {params.pre} -nt {threads} -m {params.model} "
           "-nm {params.boot_iter} -bb {params.bootstrap} > {log}"

############################
# CONSENSUS PART OF PIPELINE
############################
rule fastp:
    input:
        fw=lambda wildcards: READ_PAIRS[wildcards.basename][0],
        rv=lambda wildcards: READ_PAIRS[wildcards.basename][1],
    output:
        fw=temp('trim/{basename}/fw.fq'),
        rv=temp('trim/{basename}/rv.fq'),
        html='trim/{basename}/report.html',
        json='trim/{basename}/report.json'
    log: "log/fastp/{basename}.log"
    threads: 2
    shell:
        'fastp -i {input.fw} -I {input.rv} '
        '-o {output.fw} -O {output.rv} --html {output.html} --json {output.json} '
        '--disable_adapter_trimming --trim_poly_g --cut_tail --low_complexity_filter '
        '--complexity_threshold 50 --thread {threads} 2> {log}'

rule zip_trim:
    input:
        fw="trim/{basename}/fw.fq",
        rv="trim/{basename}/rv.fq"
    output:
        fw="trim/{basename}/fw.fq.gz",
        rv="trim/{basename}/rv.fq.gz"
    shell: "gzip -k {input}"

# Do this to get the best template, -Sparse option is designed for this.
rule initial_kma_map:
    input:
        fw=rules.fastp.output.fw,
        rv=rules.fastp.output.rv,
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

# Possibly add -ref_fsa to disallow gaps and -dense to disallow insertions
# Align to only the best template found by initial mapping round
rule kma_map:
    input:
        fw=rules.fastp.output.fw,
        rv=rules.fastp.output.rv,
        index=rules.index_ref.output,
        spa=rules.initial_kma_map.output
    output:
        res="aln/{basename}/{segment,[A-Z0-9]+}.res",
        fsa="aln/{basename}/{segment,[A-Z0-9]+}.fsa"
    params:
        db=REFOUTDIR + "/{segment}", # same as index_reffile param
        outbase="aln/{basename}/{segment}",
        subject=lambda wc: tools.get_best_subject(f"aln/{wc.basename}/{wc.segment}.spa")[1]
    log: "log/aln/{basename}_{segment}.log"
    threads: 2
    # This is a run command, because the comamdn cannot be evaluated until
    # the initial_kma_map, so printing it in a shell command will raise an error
    # If gapopen is default, it will lead to totally absurd alignments. Perhaps
    # setting open/ext/m/mm be -5/-1/1/-3 will lead to better results still?
    run:
        if params.subject is None:
            tools.touch_files(output)
        else:
            shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -gapopen -5 -nf -Mt1 " + str(params.subject) + " 2> {log}")

# In our current lab setup, we use primers to amplify our influeza segments. But these do not have the
# proper sequence. We do this before the final mapping step in order to get well-defined
# start and stop of the sequenced part for the last round of mapping.
rule remove_primers:
    input:
        con=rules.kma_map.output.fsa,
        primers=f"{SNAKEDIR}/ref/primers.fna"
    output: "seqs/{basename}/{segment}.trimmed.fna"
    log: "log/consensus/remove_primers_{basename}_{segment}.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/trim_consensus.jl",
        minmatches=6
    run:
        shell("{params.juliacmd} {params.scriptpath} {input.primers} {input.con} {output} {params.minmatches} > {log}")

# We now re-map to the created consensus sequence in order to accurately
# estimate depths and coverage, and get a more reliable assembly seq.
rule index_consensus:
    input: rules.remove_primers.output
    output:
        comp="seqs/{basename}/{segment,[a-zA-Z0-9]+}.comp.b",
        name="seqs/{basename}/{segment,[a-zA-Z0-9]+}.name",
        length="seqs/{basename}/{segment,[a-zA-Z0-9]+}.length.b",
        seq="seqs/{basename}/{segment,[a-zA-Z0-9]+}.seq.b"
    params:
        outpath="seqs/{basename}/{segment}"
    log: "log/seqs/kma_index_{basename}_{segment}.log"
    run:
        rec = tools.try_first_rec(input[0])
        if rec is None or len(rec) < 50:
            tools.touch_files(output)
        else:
            shell("kma index -i {input} -o {params.outpath} 2> {log}")

# And now we KMA map top that index again
rule final_kma_map:
    input: 
        fw=rules.fastp.output.fw,
        rv=rules.fastp.output.rv,
        index=rules.index_consensus.output,
        # We add this here so the untrimmed ones are not deleted too early.
        zipped=rules.zip_trim.output
    output:
        mat="seqs/{basename}/{segment,[A-Z0-9]+}.mat.gz",
        res="seqs/{basename}/{segment,[A-Z0-9]+}.res",
        fsa="seqs/{basename}/{segment,[A-Z0-9]+}.fsa"
    params:
        db="seqs/{basename}/{segment}",
        outbase="seqs/{basename}/{segment}",
        lengthfile="seqs/{basename}/{segment,[a-zA-Z0-9]+}.length.b"
    log: "log/seqs/kma_map_{basename}_{segment}.log"
    threads: 2
    run:
        nbytes = len(open(params.lengthfile, "rb").read())
        if nbytes < 5:
            tools.touch_files(output)
        else:
            shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -gapopen -5 -nf -matrix 2> {log}")

rule move_consensus:
    input: "seqs/{basename}/{segment}.fsa"
    output: "consensus/{basename}/{segment}.fna"
    run:
        with open(input[0], "rb") as file:
            try:
                consensus = next(tools.byte_iterfasta(file))
            except StopIteration:
                consensus = tools.FastaEntry("", bytearray())

        consensus.header = "{}_{}".format(wildcards.basename, wildcards.segment)
        with open(output[0], "w") as file:
            print(consensus.format(), file=file)

rule create_report:
    input:
        con=expand("consensus/{{basename}}/{segment}.fna", segment=SEGMENTS),
        mat=expand("seqs/{{basename}}/{segment}.mat.gz", segment=SEGMENTS),
    output:
        acc="consensus/{basename}/accepted.txt",
        rep="consensus/{basename}/report.txt",
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/report.jl"
    log: "log/report/{basename}"
    run:
        conpaths = ",".join(input.con)
        matpaths = ",".join(input.mat)
        shell(f"{{params.juliacmd}} {{params.scriptpath}} {conpaths} {matpaths} {{output}} > {{log}}")

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
    input: expand("seqs/{{basename}}/{segment}.mat.gz", segment=SEGMENTS)
    output: "depths/{basename}.pdf"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/covplot.jl"
    shell: "{params.juliacmd} {params.scriptpath} {output} {input}"

############################
# IQTREE PART OF PIPELINE
############################

def calc_subtype(wildcards):
    with open(f"aln/{wildcards.basename}/{wildcards.segment}.spa") as file:
        next(file) # skip header
        return next(file).split('\t')[0].split('|')[1]

rule cat_ref_consensus:
    input:
        consensus="consensus/{basename}/{segment}.fna",
        refaln=expand(REFOUTDIR + "/{{segment}}_{subtype}.aln.trimmed.fna", subtype=SUBTYPES)
    output: "phylogeny/{basename}/{segment}.cat.aln.trimmed.fna"
    params: lambda wc: REFOUTDIR + f"/{wc.segment}_{calc_subtype(wc)}.aln.trimmed.fna"
    log: "log/mafft/phylogeny/{basename}_{segment}.log"
    run: shell("mafft --add {input.consensus} {params} > {output} 2> {log}")

def iqtree_guide_tree(wildcards):
    subtype = calc_subtype(wildcards)
    return REFOUTDIR + f"/{wildcards.segment}_{subtype}.treefile"

rule iqtree:
    input:
        consensus=rules.cat_ref_consensus.output,
        guidetrees=expand(REFOUTDIR + "/{{segment}}_{subtype}.treefile", subtype=SUBTYPES)
    output: "phylogeny/{basename}/{segment,[A-Z0-9]+}.treefile"
    params:
        pre="phylogeny/{basename}/{segment}",
        bootstrap=1000,
        boot_iter=2500,
        model="HKY+G2",
        guide=iqtree_guide_tree
    threads: 2
    log: "log/iqtree/phylogeny/{basename}_{segment}.log"
    shell: "iqtree -s {input[0]} -pre {params.pre} -nt {threads} -m {params.model} "
           "-nm {params.boot_iter} -bb {params.bootstrap} -g {params.guide} > {log}"

"""
rule cat
_orfs:
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