import sys
import os

######################################################
# Define directories
######################################################
PROJECT_DIR = os.path.dirname(os.path.dirname(workflow.snakefile))
REF_DIR = os.path.join(PROJECT_DIR, "ref")
REF_SEQ_DIR = os.path.join(REF_DIR, "seqs")
REF_OUT_DIR = os.path.join(PROJECT_DIR, "refout")

# We have to create this directories either in a rule or outside the DAG to not
# mess up the DAG.
if not os.path.isdir(REF_OUT_DIR):
    os.mkdir(REF_OUT_DIR)

SCRIPT_DIR = os.path.join(PROJECT_DIR, "scripts")

sys.path.append(SCRIPT_DIR)
import tools

SYSIMG_PATH = os.path.join(SCRIPT_DIR, "sysimg", "sysimg.so")
JULIA_COMMAND = f"julia --startup-file=no --project={PROJECT_DIR} -J {SYSIMG_PATH}"

######################################################
# Read config
######################################################
# ------ Readdir -------
if "readdir" not in config:
    raise KeyError("You must supply read path: '--config readdir=path/to/reads'")

READ_DIR = config["readdir"]

# For some reason it all fucks up if the read directory is a child of the running
# directory, so we check that here.
_dir = os.path.abspath(READ_DIR)
while _dir != os.path.dirname(_dir):
    if _dir == os.getcwd():
        raise ValueError("Error: Read path cannot be child directory of running directory")
    _dir = os.path.dirname(_dir)

# ------ Platform ------
if "platform" not in config or config["platform"] not in ["illumina", "nanopore"]:
    raise KeyError("You must supply platform: '--config platform=[\"illumina/nanopore\"]")

IS_NANOPORE = config["platform"] == "nanopore"
IS_ILLUMINA = config["platform"] == "illumina"
assert IS_NANOPORE ^ IS_ILLUMINA

# ------ Pore ------
if IS_NANOPORE:
    if "pore" not in config:
        raise KeyError("On nanopore platform, you must supply pore: --config pore=9/10")
    if str(config["pore"]) not in ["9", "10"]:
        raise ValueError(f"Pore must be 9 or 10, not '{config['pore']}'")
    PORE = str(config["pore"])

#################################
# Rules
#################################
rule index_ref:
    input: REF_SEQ_DIR + "/{segment}.fna"
    output:
        comp=REF_OUT_DIR + "/{segment,[a-zA-Z0-9]+}.comp.b",
        name=REF_OUT_DIR + "/{segment,[a-zA-Z0-9]+}.name",
        length=REF_OUT_DIR + "/{segment,[a-zA-Z0-9]+}.length.b",
        seq=REF_OUT_DIR + "/{segment,[a-zA-Z0-9]+}.seq.b"
    params:
        outpath=REF_OUT_DIR + "/{segment}",
    log: "log/kma_ref/{segment}.log"
    
    shell: "kma index -k 12 -Sparse - -i {input} -o {params.outpath} 2> {log}"
