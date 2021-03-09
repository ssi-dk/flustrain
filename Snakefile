import sys
import os
import itertools as it

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools

# We only need this because Julia 1.6 segfaults with PkgCompiler at the moment
SYSIMG_PATH = os.path.join(SNAKEDIR, "scripts", "sysimg", "sysimg.so")
JULIA_COMMAND = f"julia --startup-file=no --project={SNAKEDIR} -J {SYSIMG_PATH}"

######################################################
# GLOBAL CONSTANTS
######################################################
if "readdir" not in config:
    raise KeyError("You must supply absolute read path: '--config readdir=/path/to/reads'")

READDIR = config["readdir"]

if "refset" not in config:
    raise KeyError("You must supply name of reference set: '--config refset=swine'")

if config["refset"] not in ["human", "swine", "avian"]:
    raise KeyError("refset must be 'human', 'swine' or 'avian'")

# For some reason it all fucks up if the read directory is a child of the running
# directory, so we check that here.
_dir = os.path.abspath(READDIR)
while _dir != os.path.dirname(_dir):
    if _dir == os.getcwd():
        raise ValueError("Error: Read path cannot be child directory of running directory")
    _dir = os.path.dirname(_dir)

READ_PAIRS = tools.get_read_pairs(READDIR)
BASENAMES = sorted(READ_PAIRS.keys())

REFSEQDIR = os.path.join(SNAKEDIR, "ref", "seqs")
REFTREEDIR = os.path.join(SNAKEDIR, "ref", "trees", config["refset"])
for path in REFSEQDIR, REFTREEDIR:
    if not os.path.isdir(path):
        raise NotADirectoryError(f"Could not find reference directory: '{path}'")

REFOUTDIR = os.path.join(SNAKEDIR, "refout")
REFOUTTREEDIR = os.path.join(REFOUTDIR, config["refset"])

# We have to create these directories either in a rule or outside the DAG to not
# mess up the DAG.
for path in [REFOUTDIR, REFOUTTREEDIR]:
    if not os.path.isdir(path):
        os.mkdir(path)

SEGMENTS = sorted([fn.rpartition('.')[0] for fn in os.listdir(REFSEQDIR) if fn.endswith(".fna")])

TREE_SEGMENTS = [fn.partition(".")[0] for fn in os.listdir(REFTREEDIR)]

if not set(TREE_SEGMENTS).issubset(set(SEGMENTS)):
    raise ValueError(f"Tree segment name in {REFTREEDIR} must be subset of segment names.")

######################################################
# Start of pipeline
######################################################
ruleorder: align_references > translate > trim_alignment > mafft > second_kma_map > first_kma_map > gzip

# We add the checkpoiint output here just to make sure the checkpoint is included
# in the DAG


def done_input(wildcards):
    # Add report and the commit
    inputs = ["report.txt"]

    # By adding this checkpoint, we can access stuff not existing yet
    checkpoints.create_report.get()

    for basename in BASENAMES:
        # Add trimmed FASTQ
        for direction in ["fw", "rv"]:
            inputs.append(f"trim/{basename}/{direction}.fq.gz")

        # Add depth plot
        inputs.append(f"depths/{basename}.pdf")
        inputs.append(f"consensus/{basename}/consensus.fna")

        # Add phylogeny - only for quality-controlled tree segments
        with open(f"consensus/{basename}/curated.fna") as file:
            for line in file:
                if line.startswith(">"):
                    segment = line.strip().rpartition("_")[-1]
                    assert segment in SEGMENTS
                    if segment in TREE_SEGMENTS:
                        inputs.append(f"phylogeny/{basename}/{segment}.treefile")

    # TODO: Add mutation scanning here.

    return inputs

rule all:
    input: done_input
    output: "commit.txt"
    params: SNAKEDIR
    shell: "git -C {params} rev-parse --short HEAD > {output} && cp {params}/copy_readme.md README.md"

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

"""
rule extract_orf:
    input: "{dir}/{base}.fna"
    output: "{dir}/{base}.orf.fna"
    run:
        with open(input[0], "rb") as file, open(output[0], "w") as outfile:
            for entry in tools.byte_iterfasta(file):
                pos, orf = tools.find_orf(entry)
                print(orf.format(), file=outfile)
"""

#################################
# REFERENCE-ONLY PART OF PIPELINE
#################################
rule index_ref:
    input: REFSEQDIR + "/{segment}.fna"
    output:
        comp=REFOUTDIR + "/{segment,[a-zA-Z0-9]+}.comp.b",
        name=REFOUTDIR + "/{segment,[a-zA-Z0-9]+}.name",
        length=REFOUTDIR + "/{segment,[a-zA-Z0-9]+}.length.b",
        seq=REFOUTDIR + "/{segment,[a-zA-Z0-9]+}.seq.b"
    params:
        outpath=REFOUTDIR + "/{segment}",
    log: "log/kma_ref/{segment}.log"
    
    shell: "kma index -k 12 -Sparse - -i {input} -o {params.outpath} 2> {log}"

rule align_references:
    input: REFTREEDIR + "/{segment}.fna"
    output: REFOUTTREEDIR + "/{segment}.aln.fna"
    threads: 4
    log: "log/align_ref/{segment}.log"
    shell: "mafft --thread {threads} --auto {input} > {output} 2> {log}"

rule reference_iqtree:
    input: REFOUTTREEDIR + "/{segment}.aln.trimmed.fna"
    output: REFOUTTREEDIR + "/{segment}.treefile"
    params:
        pre=REFOUTTREEDIR + "/{segment}",
        model="HKY+G2" # just pick a model to be consistent
    threads: 2
    log: "log/iqtree/{segment}.log"
    shell: "iqtree -s {input} -pre {params.pre} -nt {threads} -m {params.model} > {log}"

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

rule gzip:
    input: "{base}"
    output: "{base}.gz"
    shell: "gzip -k {input}"

# Do this to get the best template, -Sparse option is designed for this.
rule map_best_template:
    input:
        fw=rules.fastp.output.fw,
        rv=rules.fastp.output.rv,
        index=rules.index_ref.output
    output: "aln/{basename}/{segment,[A-Z0-9]+}.spa"
    params:
        db=REFOUTDIR + "/{segment}", # same as index_reffile param
        outbase="aln/{basename}/{segment}"
    threads: 1
    log: "log/aln/{basename}_{segment}.initial.log"
    shell:
        # Sort by template coverate (-ss c), otherwise "best" template will
        # be driven by erroneous sequencing kmers
        "kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
        "-t {threads} -Sparse -ss c 2> {log}"

rule collect_best_templates:
    input: expand("aln/{basename}/{segment}.spa", segment=SEGMENTS, basename=BASENAMES)
    output: temp(expand("aln/{basename}/cat.fna", basename=BASENAMES))
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/gather_spa.jl",
        refpath=REFSEQDIR
    shell: "julia {params.scriptpath} aln {params.refpath}"

rule first_kma_index:
    input: "aln/{basename}/cat.fna"
    output:
        comp=temp("aln/{basename}/cat.comp.b"),
        name=temp("aln/{basename}/cat.name"),
        length=temp("aln/{basename}/cat.length.b"),
        seq=temp("aln/{basename}/cat.seq.b")
    params:
        t_db="aln/{basename}/cat"
    log: "log/aln/kma1_index_{basename}.log"
    # The pipeline is very sensitive to the value of k here.
    # Too low means the mapping is excruciatingly slow,
    # too high results in poor mapping quality.
    # k = 10 might be a little on the low side.
    shell: "kma index -k 10 -i {input} -o {params.t_db} 2> {log}"

rule first_kma_map:
    input:
        fw=rules.fastp.output.fw,
        rv=rules.fastp.output.rv,
        index=rules.first_kma_index.output,
    output:
        res="aln/{basename}/kma1.res",
        fsa="aln/{basename}/kma1.fsa",
        mat="aln/{basename}/kma1.mat.gz",
    params:
        db="aln/{basename}/cat",
        outbase="aln/{basename}/kma1",
    log: "log/aln/kma1_map_{basename}.log"
    threads: 2
    run:
        shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
        "-t {threads} -1t1 -gapopen -5 -nf -matrix 2> {log}")

# In our current lab setup, we use primers to amplify our influeza segments. But these do not have the
# proper sequence. We do this before the final mapping step in order to get well-defined
# start and stop of the sequenced part for the last round of mapping.

rule remove_primers:
    input:
        con=rules.first_kma_map.output.fsa,
        primers=f"{SNAKEDIR}/ref/primers.fna"
    output: temp("aln/{basename}/cat.trimmed.fna")
    log: "log/consensus/remove_primers_{basename}.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/trim_consensus.jl",
        minmatches=4,
        fuzzylen=8,
    run:
        shell("{params.juliacmd} {params.scriptpath} {input.primers} "
              "{input.con} {output} {params.minmatches} {params.fuzzylen} > {log}")

# We now re-map to the created consensus sequence in order to accurately
# estimate depths and coverage, and get a more reliable assembly seq.
rule second_kma_index:
    input: "aln/{basename}/cat.trimmed.fna"
    output:
        comp=temp("aln/{basename}/cat.trimmed.comp.b"),
        name=temp("aln/{basename}/cat.trimmed.name"),
        length=temp("aln/{basename}/cat.trimmed.length.b"),
        seq=temp("aln/{basename}/cat.trimmed.seq.b")
    params:
        t_db="aln/{basename}/cat.trimmed"
    log: "log/aln/kma2_index_{basename}.log"
    shell: "kma index -i {input} -o {params.t_db} 2> {log}"

# And now we KMA map to that index again
rule second_kma_map:
    input:
        fw=rules.fastp.output.fw,
        rv=rules.fastp.output.rv,
        index=rules.second_kma_index.output,
        keepfw="trim/{basename}/fw.fq.gz",
        keeprv="trim/{basename}/rv.fq.gz"
    output:
        res="aln/{basename}/kma2.res",
        fsa="aln/{basename}/kma2.fsa"
    params:
        db="aln/{basename}/cat.trimmed",
        outbase="aln/{basename}/kma2",
    log: "log/aln/kma2_map_{basename}.log"
    threads: 2
    run:
        shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
        "-t {threads} -1t1 -gapopen -5 -nf 2> {log}")

checkpoint create_report:
    input:
        matrix=expand("aln/{basename}/kma1.mat.gz", basename=BASENAMES),
        assembly=expand("aln/{basename}/kma2.fsa", basename=BASENAMES),
        res=expand("aln/{basename}/kma2.res", basename=BASENAMES)
    output:
        consensus=expand("consensus/{basename}/{type}.{nuc}",
            basename=BASENAMES, type=["consensus", "curated"], nuc=["fna", "faa"]),
        report="report.txt",
        depths=expand("depths/{basename}.pdf", basename=BASENAMES)
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/report.jl",
        refdir=REFSEQDIR
    log: "log/report.txt"
    threads: workflow.cores
    run:
        shell(f"julia -t {threads} {params.scriptpath} consensus report.txt "
               "depths aln kma2.fsa kma2.res kma1.mat.gz {params.refdir} > {log}")

############################
# IQTREE PART OF PIPELINE
############################
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

rule cat_ref_consensus:
    input:
        consensus="phylogeny/{basename}/{segment}.fna",
        refaln=REFOUTTREEDIR + "/{segment}.aln.trimmed.fna"
    output: temp("phylogeny/{basename}/{segment}.cat.aln.trimmed.fna")
    log: "log/mafft/phylogeny/{basename}_{segment}.log"
    run: shell("mafft --add {input.consensus} {input.refaln} > {output} 2> {log}")

rule iqtree:
    input:
        consensus=rules.cat_ref_consensus.output,
        guidetree=rules.reference_iqtree.output
    output: "phylogeny/{basename}/{segment,[A-Z0-9]+}.treefile"
    params:
        pre="phylogeny/{basename}/{segment}",
        model="HKY+G2",
        guide=REFOUTTREEDIR + "/{segment}.treefile"
    threads: 2
    log: "log/iqtree/phylogeny/{basename}_{segment}.log"
    shell: "iqtree -s {input[0]} -pre {params.pre} -nt {threads} -m {params.model} "
           "-g {params.guide} > {log}"

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
