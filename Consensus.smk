import sys
import os

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools

# We only need this because Julia 1.6 segfaults with PkgCompiler at the moment
SYSIMG_PATH = os.path.join(SNAKEDIR, "scripts", "sysimg", "sysimg.so")

# TODO: I've removed -J {SYSIMG_PATH} here while I rewrite the workflow. Can be re-added later.
JULIA_COMMAND = f"julia --startup-file=no --project={SNAKEDIR}"

######################################################
# GLOBAL CONSTANTS
######################################################
# ------ Readdir -------
if "readdir" not in config:
    raise KeyError("You must supply read path: '--config readdir=path/to/reads'")

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
assert IS_NANOPORE ^ IS_ILLUMINA # only one must be true at a time

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

# TODO: After report.jl refactor, check all outputs are present here.
def done_input(wildcards):
    # Add report and the commit
    inputs = ["report.txt"]

    # Add trimmed FASTQ
    for basename in BASENAMES:
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
    log: "tmp/log/kma_ref/{segment}.log"
    
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
            fw=temp('tmp/trim/{basename}/fw.fq'),
            rv=temp('tmp/trim/{basename}/rv.fq'),
            html='tmp/trim/{basename}/report.html',
            json='tmp/trim/{basename}/report.json'
        log: "tmp/log/fastp/{basename}.log"
        threads: 2
        shell:
            'fastp -i {input.fw} -I {input.rv} '
            '-o {output.fw} -O {output.rv} --html {output.html} --json {output.json} '
            '--disable_adapter_trimming --trim_poly_g --cut_tail --cut_front --low_complexity_filter '
            '--complexity_threshold 50 --thread {threads} 2> {log}'

    rule map_best_template:
        input:
            fw=rules.fastp.output.fw,
            rv=rules.fastp.output.rv,
            index=rules.index_ref.output
        output: "tmp/aln/{basename}/{segment,[A-Z0-9]+}.spa"
        params:
            db=REFOUTDIR + "/{segment}", # same as index_reffile param
            outbase="tmp/aln/{basename}/{segment}"
        threads: 1
        log: "tmp/log/aln/{basename}_{segment}.initial.log"
        shell:
            # Here, not sure if I should sort by template cov (-ss c) or not.
            # Pro: We mostly care about having a fully-covered template, not a partially
            # covered with high depth at the covered areas:
            # Con: Majority vote will win away, so it'll just fuck up if we pick a
            # uniformly, but low covered reference anyway
            "kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
            "-ss c -t {threads} -Sparse 2> {log}"

elif IS_NANOPORE:
    rule fastp:
        input: lambda wildcards: READS[wildcards.basename]
        output:
            reads=temp('tmp/trim/{basename}/reads.fq'),
            html='tmp/trim/{basename}/report.html',
            json='tmp/trim/{basename}/report.json'
        log: "tmp/log/fastp/{basename}.log"
        threads: 2
        shell:
            'fastp -i {input} -o {output.reads} --html {output.html} '
            '--json {output.json} --disable_adapter_trimming  --disable_trim_poly_g '
            '--cut_window_size 10 --cut_mean_quality 10 --cut_tail --cut_front --low_complexity_filter  '
            '--complexity_threshold 50 --length_limit 2400 --length_required 100 '
            '--average_qual 12 --thread {threads} 2> {log}'

    rule map_best_template:
        input:
            reads=rules.fastp.output.reads,
            index=rules.index_ref.output
        output: "tmp/aln/{basename}/{segment,[A-Z0-9]+}.spa"
        params:
            db=REFOUTDIR + "/{segment}", # same as index_reffile param
            outbase="tmp/aln/{basename}/{segment}"
        threads: 1
        log: "tmp/log/aln/{basename}_{segment}.initial.log"
        shell:
            # See above comment in rule with same name
            "kma -i {input.reads} -o {params.outbase} -t_db {params.db} "
            "-ss c -t {threads} -Sparse 2> {log}"

### Both platforms
rule collect_best_templates:
    input: expand("tmp/aln/{basename}/{segment}.spa", segment=SEGMENTS, basename=BASENAMES)
    output: expand("tmp/aln/{basename}/cat.fna", basename=BASENAMES)
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/gather_spa.jl",
        refpath=REFSEQDIR
    shell: "{params.juliacmd} {params.scriptpath} tmp/aln {params.refpath}"

rule first_kma_index:
    input: "tmp/aln/{basename}/cat.fna"
    output:
        comp=temp("tmp/aln/{basename}/cat.comp.b"),
        name=temp("tmp/aln/{basename}/cat.name"),
        length=temp("tmp/aln/{basename}/cat.length.b"),
        seq=temp("tmp/aln/{basename}/cat.seq.b")
    params:
        t_db="tmp/aln/{basename}/cat"
    log: "tmp/log/aln/kma1_index_{basename}.log"
    # The pipeline is very sensitive to the value of k here.
    # Too low means the mapping is excruciatingly slow,
    # too high results in poor mapping quality.
    # k = 10 might be a little on the low side.
    shell: "kma index -k 10 -i {input} -o {params.t_db} 2> {log}"

if IS_ILLUMINA:
    rule first_kma_map:
        input:
            fw=rules.fastp.output.fw,
            rv=rules.fastp.output.rv,
            index=rules.first_kma_index.output,
        output:
            res="tmp/aln/{basename}/kma1.res",
            fsa="tmp/aln/{basename}/kma1.fsa",
            mat="tmp/aln/{basename}/kma1.mat.gz",
        params:
            db="tmp/aln/{basename}/cat",
            outbase="tmp/aln/{basename}/kma1",
        log: "tmp/log/tmp/aln/kma1_map_{basename}.log"
        threads: 2
        run:
            shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -1t1 -gapopen -5 -nf -matrix 2> {log}")

elif IS_NANOPORE:
    rule first_kma_map:
        input:
            reads=rules.fastp.output.reads,
            index=rules.first_kma_index.output,
        output:
            res="tmp/aln/{basename}/kma1.res",
            fsa="tmp/aln/{basename}/kma1.fsa",
            mat="tmp/aln/{basename}/kma1.mat.gz",
        params:
            db="tmp/aln/{basename}/cat",
            outbase="tmp/aln/{basename}/kma1",
        log: "tmp/log/aln/kma1_map_{basename}.log"
        threads: 2
        run:
            shell("kma -i {input.reads} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -1t1 -bcNano -nf -matrix 2> {log}")

# Both platforms
rule remove_primers:
    input:
        con=rules.first_kma_map.output.fsa,
        primers=f"{SNAKEDIR}/ref/primers.fna"
    output: temp("tmp/aln/{basename}/cat.trimmed.fna")
    log: "tmp/log/consensus/remove_primers_{basename}.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/trim_consensus.jl",
        minmatches=4,
        fuzzylen=8,
    shell:
        "{params.juliacmd} {params.scriptpath} {input.primers} "
        "{input.con} {output} {params.minmatches} {params.fuzzylen} > {log}"

if IS_ILLUMINA:
    # We now re-map to the created consensus sequence in order to accurately
    # estimate depths and coverage, and get a more reliable assembly seq.
    rule second_kma_index:
        input: rules.remove_primers.output
        output:
            comp=temp("tmp/aln/{basename}/cat.trimmed.comp.b"),
            name=temp("tmp/aln/{basename}/cat.trimmed.name"),
            length=temp("tmp/aln/{basename}/cat.trimmed.length.b"),
            seq=temp("tmp/aln/{basename}/cat.trimmed.seq.b")
        params:
            t_db="tmp/aln/{basename}/cat.trimmed"
        log: "tmp/log/aln/kma2_index_{basename}.log"
        shell: "kma index -i {input} -o {params.t_db} 2> {log}"

    # And now we KMA map to that index again
    rule second_kma_map:
        input:
            fw=rules.fastp.output.fw,
            rv=rules.fastp.output.rv,
            index=rules.second_kma_index.output,
        output:
            res="tmp/aln/{basename}/kma2.res",
            fsa="tmp/aln/{basename}/kma2.fsa",
            mat="tmp/aln/{basename}/kma2.mat.gz"
        params:
            db="tmp/aln/{basename}/cat.trimmed",
            outbase="tmp/aln/{basename}/kma2",
        log: "tmp/log/aln/kma2_map_{basename}.log"
        threads: 2
        run:
            shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -1t1 -gapopen -5 -nf -matrix 2> {log}")

    rule create_report:
        input:
            matrix=expand("tmp/aln/{basename}/kma1.mat.gz", basename=BASENAMES),
            assembly=expand("tmp/aln/{basename}/kma2.fsa", basename=BASENAMES),
            res=expand("tmp/aln/{basename}/kma2.res", basename=BASENAMES)
        output:
            consensus=expand("consensus/{basename}/{type}.{nuc}",
                basename=BASENAMES, type=["consensus", "curated"], nuc=["fna", "faa"]
            ),
            report="report.txt",
            depths=expand("depths/{basename}.pdf", basename=BASENAMES)
        params:
            juliacmd=JULIA_COMMAND,
            scriptpath=f"{SNAKEDIR}/scripts/report.jl",
            refdir=REFSEQDIR
        log: "tmp/log/report.txt"
        threads: workflow.cores
        run:
            shell(f"{params.juliacmd} -t {threads} {params.scriptpath} illumina . {params.refdir} > {log}")

elif IS_NANOPORE:
    rule medaka:
        input: 
            reads=rules.fastp.output.reads,
            draft=rules.remove_primers.output
        output: directory("tmp/aln/{basename}/medaka")
        log: "tmp/log/aln/medaka_{basename}.log"
        threads: 2
        params:
            model=lambda wc: "r941_min_high_g360" if PORE == 9 else "r103_min_high_g360"
        shell:
            "medaka_consensus -i {input.reads} -d {input.draft} -o {output} -t {threads} "
            "-m {params.model} 2> {log}"

    rule clean_medaka:
        input: rules.medaka.output
        output: "tmp/aln/{basename}/moved.txt"
        shell: "rm {input}/*.bam {input}/*.bai {input}/*.hdf && touch {output}"

    rule create_report:
        input:
            assembly=expand("tmp/aln/{basename}/moved.txt", basename=BASENAMES),
        output:
            consensus=expand("consensus/{basename}/{type}.{nuc}",
                basename=BASENAMES, type=["consensus", "curated"], nuc=["fna", "faa"]
            ),
            report="report.txt",
        params:
            juliacmd=JULIA_COMMAND,
            scriptpath=f"{SNAKEDIR}/scripts/InfluenzaReport/src/InfluenzaReport.jl",
            refdir=REFSEQDIR
        log: "tmp/log/report.txt"
        threads: workflow.cores
        run:
            shell(f"{params.juliacmd} -t {threads} {params.scriptpath} nanopore . {params.refdir} > {log}")
