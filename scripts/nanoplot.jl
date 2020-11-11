using Plots
using FASTX
using CodecZlib
using Comonicon

function maybe_decompress(path)
    io = open(path)
    sig = read(io, 2)
    seekstart(io)
    if length(sig) < 2 || sig[1:2] != [0x1f, 0x8b]
        return io
    else
        return GzipDecompressorStream(io)
    end
end

# This converts from PHRED33 to a probability
qual(x::UInt8) = 10^-((x - UInt8(33)) / 10)

function meanqual(rec::FASTQ.Record)
    sumqual = zero(Float64)
    data = rec.data
    # We need to average by averaging the float probabilities,
    # then converting back to phred
    @inbounds for i in rec.quality
        sumqual += qual(data[i])
    end
    meanqual = sumqual / length(rec.quality)
    # Convert back to PHRED
    phredmean = -10 * log10(meanqual)

    return length(rec.quality), phredmean
end

function getstats(fastq::IO)
    reader = FASTQ.Reader(fastq)
    record = FASTQ.Record()
    meanquals = Float32[]
    lengths = UInt32[]
    while !eof(reader)
        read!(reader, record)
        len, qual = meanqual(record)
        push!(lengths, len)
        push!(meanquals, qual)
    end
    return lengths, meanquals
end

function cumulative(lengths)
    totalbases = sum(lengths)
    bases = totalbases
    slens = sort(lengths)
    ys = Float32[bases]
    i, minlength, len = 1, 0, 0
    while true
        while len < minlength
            i += 1
            i > length(slens) && break
            len = slens[i]
            bases -= len
        end
        push!(ys, bases)
        i > length(slens) && break
        # Once we drop below 1%, we don't care for the rest. 
        bases < 0.01 * totalbases && break
        minlength += 10
    end
    return ys, range(0, step=10, length=length(ys))
end

hist(y, xlabel) = histogram(y, xlabel=xlabel, legend=nothing)

"""
Plot quality/length statistics for a PHRED33 FASTQ file.

# Arguments

- `inpath`: Input path to (gzipped or plain) FASTQ file
- `outdir`: Output directory to create
"""
@main function main(inpath::String, outdir::String)
    mkdir(outdir)
    io = maybe_decompress(inpath)
    lengths, meanquals = getstats(io)
    close(io)

    # Histogram of qualities
    plt = hist(meanquals, "Mean PHRED score")
    savefig(plt, joinpath(outdir, "quals.png"))

    # Read lengths
    plt = hist(log10.(lengths), "Log10 read length")
    savefig(plt, joinpath(outdir, "lengths.log.png"))

    # Weighted read lengths
    weighted_lens = lengths .* meanquals
    plt = hist(log10.(weighted_lens), "Log10 read length (weighted)")
    savefig(plt, joinpath(outdir, "lengths.weighted.log.png"))
    weighted_lens = nothing

    # Cumulative
    ys, xs = cumulative(lengths)
    plt = plot(xs, ys, xlabel="Minimum read size", ylabel="Cumulative bases")
    savefig(plt, joinpath(outdir, "lengths.cumulative.png"))
end