from snakemake.utils import min_version
min_version("6.0")

module consflow:
    snakefile: "Consensus.smk"

use rule * from consflow as consflow_*

module phyloflow:
    snakefile: "Phylogeny.smk"

use rule * from phyloflow as phyloflow_*

print("foo")