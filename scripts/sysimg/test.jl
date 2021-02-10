REFDIR = "/Users/jakobnissen/Documents/ssi/projects/flupipe/ref/swine"

module Trim
    include("../trim_consensus.jl")
end

Trim.trim_consensus("primers.fna", "test.fna", "out.fna", 4)

module Report
    include("../report.jl")
end

mkdir("out")
mkdir("out/depths")
mkdir("out/foo")
Report.main("out", "out/report.txt", "out/depths", "aln", "kma2.fsa", "kma2.mat.gz", REFDIR)

#=
module Nanofilt
    include("../nanofilt.jl")
end

Nanofilt._nanofilt("test.fastq.gz", "test.filt.gz", 10, 50, 20)

module Nanoplot
    include("../nanoplot.jl")
end

Nanoplot._nanoplot("test.fastq.gz", "out")
=#
