include: "consensus_config.snakefile"

# Constants
READS = tools.get_read_pairs(READ_DIR)
BASENAMES = sorted(READS.keys())
SEGMENTS = sorted([fn.rpartition('.')[0] for fn in os.listdir(REF_SEQ_DIR) if fn.endswith(".fna")])

# TODO: After report.jl refactor, check all outputs are present here.
def done_input(wildcards):
    # Add report and the commit
    inputs = ["report.txt"]

    for basename in BASENAMES:
        # Add trimmed FASTQ
        for direction in ["fw", "rv"]:
            inputs.append(f"trim/{basename}/{direction}.fq.gz")

    return inputs

ruleorder: second_kma_map > first_kma_map > gzip

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
        db=REF_OUT_DIR + "/{segment}", # same as index_reffile param
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
        scriptpath=f"{SCRIPT_DIR}/gather_spa.jl",
        refpath=REF_SEQ_DIR
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

rule remove_primers:
    input:
        con=rules.first_kma_map.output.fsa,
        primers=f"{REF_DIR}/primers.fna"
    output: temp("aln/{basename}/cat.trimmed.fna")
    log: "log/consensus/remove_primers_{basename}.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SCRIPT_DIR}/trim_consensus.jl",
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
        scriptpath=f"{SCRIPT_DIR}/report.jl",
        refdir=REF_SEQ_DIR
    log: "log/report.txt"
    threads: workflow.cores
    run:
        shell(f"julia -t {threads} {params.scriptpath} consensus report.txt "
               "depths aln kma2.fsa kma2.res kma1.mat.gz {params.refdir} > {log}")