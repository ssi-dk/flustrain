module Trim
    include("../trim_consensus.jl")
end

Trim.trim_consensus("primers.fna", "test.fna", "out.fna", 4)

module CovPlot
    include("../covplot.jl")
end

CovPlot.main("out.pdf", ["NS.mat.gz"])

module Report
    include("../report.jl")
end

Report.main(["NS.fsa"], ["NS.mat.gz"], "out.acc", "out.report")

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
