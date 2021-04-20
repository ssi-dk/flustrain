include: "consensus_config.snakefile"

# Constants
READS = tools.get_nanopore_reads(READ_DIR)
BASENAMES = sorted(READS.keys())
SEGMENTS = sorted([fn.rpartition('.')[0] for fn in os.listdir(REF_SEQ_DIR) if fn.endswith(".fna")])

# TODO: After report.jl refactor, check all outputs are present here.
def done_input(wildcards):
    # Add report and the commit
    #inputs = ["report.txt"]
    inputs = []

    for basename in BASENAMES:
        # Add trimmed FASTQ
        inputs.append(f"trim/{basename}/reads.fq.gz")
        inputs.append(f"aln/{basename}/medaka")

    return inputs

rule all:
    input: done_input
    output: "commit.txt"
    params: PROJECT_DIR
    shell: "git -C {params} rev-parse --short HEAD > {output} && cp {params}/copy_readme.md README.md"

rule gzip:
    input: "{base}"
    output: "{base}.gz"
    shell: "gzip -k {input}"

rule fastp:
    input: lambda wildcards: READS[wildcards.basename]
    output:
        reads=temp('trim/{basename}/reads.fq'),
        html='trim/{basename}/report.html',
        json='trim/{basename}/report.json'
    log: "log/fastp/{basename}.log"
    threads: 2
    shell:
        'fastp -i {input} -o {output.reads} --html {output.html} '
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
        db=REF_OUT_DIR + "/{segment}", # same as index_reffile param
        outbase="aln/{basename}/{segment}"
    threads: 1
    log: "log/aln/{basename}_{segment}.initial.log"
    shell:
        # Sort by template coverate (-ss c), otherwise "best" template will
        # be driven by erroneous sequencing kmers
        "kma -i {input.reads} -o {params.outbase} -t_db {params.db} "
        "-t {threads} -Sparse -ss c 2> {log}"

rule collect_best_templates:
    input: expand("aln/{basename}/{segment}.spa", segment=SEGMENTS, basename=BASENAMES)
    output: temp(expand("aln/{basename}/cat.fna", basename=BASENAMES))
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SCRIPT_DIR}/gather_spa.jl",
        refpath=REF_SEQ_DIR
    shell: "julia {params.scriptpath} aln {params.refpath}"

rule medaka:
    input: 
        reads=rules.fastp.output.reads,
        draft="aln/{basename}/cat.fna"
    output: directory("aln/{basename}/medaka")
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
        primers=f"{REF_DIR}/primers.fna"
    output: temp("aln/{basename}/consensus.trimmed.fna")
    log: "log/consensus/remove_primers_{basename}.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SCRIPT_DIR}/trim_consensus.jl",
        minmatches=4,
        fuzzylen=8,
    run:
        shell("{params.juliacmd} {params.scriptpath} {input.primers} "
            "{input.con}/consensus.fasta {output} {params.minmatches} {params.fuzzylen} > {log}")
