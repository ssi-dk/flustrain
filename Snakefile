import sys
import os

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools

# We only need this because Julia 1.6 segfaults with PkgCompiler at the moment
SYSIMG_PATH = os.path.join(SNAKEDIR, "scripts", "sysimg", "sysimg.so")
JULIA_COMMAND = f"julia --startup-file=no --project={SNAKEDIR} -J {SYSIMG_PATH}"

######################################################
# GLOBAL CONSTANTS
######################################################
# ------ Readdir -------
if "readdir" not in config:
    raise KeyError("You must supply absolute read path: '--config readdir=/path/to/reads'")

READDIR = config["readdir"]

# For some reason it all fucks up if the read directory is a child of the running
# directory, so we check that here.
_dir = os.path.abspath(READDIR)
while _dir != os.path.dirname(_dir):
    if _dir == os.getcwd():
        raise ValueError("Error: Read path cannot be child directory of running directory")
    _dir = os.path.dirname(_dir)

if "platform" not in config or config["platform"] not in ["illumina", "nanopore"]:
    raise KeyError("You must supply platform: '--config platform=[\"illumina/nanopore\"]")

IS_NANOPORE = config["platform"] == "nanopore"
IS_ILLUMINA = config["platform"] == "illumina"
assert IS_NANOPORE ^ IS_ILLUMINA

if IS_NANOPORE:
    if "pore" not in config:
        raise KeyError("On nanopore platform, you must supply pore: --config pore=9/10")
    if str(config["pore"]) not in ["9", "10"]:
        raise ValueError(f"Pore must be 9 or 10, not '{config['pore']}'")
    PORE = str(config["pore"])

if IS_ILLUMINA:
    READS = tools.get_read_pairs(READDIR)
else:
    READS = tools.get_nanopore_reads(READDIR)

BASENAMES = sorted(READS.keys())

REFSEQDIR = os.path.join(SNAKEDIR, "ref", "seqs")
REFOUTDIR = os.path.join(SNAKEDIR, "refout")

# We have to create this directories either in a rule or outside the DAG to not
# mess up the DAG.
if not os.path.isdir(REFOUTDIR):
    os.mkdir(REFOUTDIR)

SEGMENTS = sorted([fn.rpartition('.')[0] for fn in os.listdir(REFSEQDIR) if fn.endswith(".fna")])

######################################################
# Start of pipeline
######################################################

# We add the checkpoiint output here just to make sure the checkpoint is included
# in the DAG

# TODO: After report.jl refactor, check all outputs are present here.
def done_input(wildcards):
    # Add report and the commit
    inputs = ["report.txt"]

    # By adding this checkpoint, we can access stuff not existing yet
    checkpoints.create_report.get()

    for basename in BASENAMES:
        # Add trimmed FASTQ
        for direction in ["fw", "rv"]:
            inputs.append(f"trim/{basename}/{direction}.fq.gz")

        inputs.append(f"consensus/{basename}/consensus.fna")

    return inputs

rule all:
    input: done_input
    output: "commit.txt"
    params: SNAKEDIR
    shell: "git -C {params} rev-parse --short HEAD > {output} && cp {params}/copy_readme.md README.md"

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

############################
# CONSENSUS PART OF PIPELINE
############################
rule gzip:
    input: "{base}"
    output: "{base}.gz"
    shell: "gzip -k {input}"

if IS_ILLUMINA:
    ruleorder: second_kma_map > first_kma_map > gzip
    rule fastp:
        input:
            fw=lambda wildcards: READS[wildcards.basename][0],
            rv=lambda wildcards: READS[wildcards.basename][1],
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

elif IS_NANOPORE:
    rule fastp:
        input: lambda wildcards: READS[wildcards.basename]
        output:
            reads=temp('trim/{basename}/reads.fq'),
            html='trim/{basename}/report.html',
            json='trim/{basename}/report.json'
        log: "log/fastp/{basename}.log"
        threads: 2
        shell:
            'fastp -i {input} -o {output.reads} -O {output.rv} --html {output.html} '
            '--json {output.json} --disable_adapter_trimming  --disable_trim_poly_g '
            '--cut_window_size 10 --cut_mean_quality 10 --low_complexity_filter  '
            '--complexity_threshold 40 --length_limit 2400 --length_required 100 '
            '--average_qual 12 --thread {threads} 2> {log}'

    rule map_best_template:
        input:
            reads=rules.fastp.output.reads,
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
            "kma -i {input.reads} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -Sparse -ss c 2> {log}"

### Both platforms
rule collect_best_templates:
    input: expand("aln/{basename}/{segment}.spa", segment=SEGMENTS, basename=BASENAMES)
    output: temp(expand("aln/{basename}/cat.fna", basename=BASENAMES))
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/gather_spa.jl",
        refpath=REFSEQDIR
    shell: "julia {params.scriptpath} aln {params.refpath}"

if IS_ILLUMINA:
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
elif IS_NANOPORE:
    rule medaka:
        input: 
            reads=rules.fastp.output.reads,
            draft="aln/{basename}/cat.fna"
        output: "aln/{basename}/medaka"
        log: "log/aln/medaka_{basename}.log"
        threads: 2
        params:
            model=lambda wc: "r941_min_high_g360" if PORE == 9 else "r103_min_high_g360"
        shell: 
            "medaka_consensus -i {input.reads} -d {input.draft} -o {output} -t {threads} "
            "-m {params.model}"

    # In our current lab setup, we use primers to amplify our influeza segments. But these do not have the
    # proper sequence. We do this before the final mapping step in order to get well-defined
    # start and stop of the sequenced part for the last round of mapping.
    rule remove_primers:
        input:
            con=rules.medaka.output,
            primers=f"{SNAKEDIR}/ref/primers.fna"
        output: temp("aln/{basename}/consensus.trimmed.fna")
        log: "log/consensus/remove_primers_{basename}.txt"
        params:
            juliacmd=JULIA_COMMAND,
            scriptpath=f"{SNAKEDIR}/scripts/trim_consensus.jl",
            minmatches=4,
            fuzzylen=8,
        run:
            shell("{params.juliacmd} {params.scriptpath} {input.primers} "
                "{input.con}/consensus.fasta {output} {params.minmatches} {params.fuzzylen} > {log}")

### Both platforms
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