using FASTX
using CodecZlib
using Comonicon

struct Options
    minlength::Int
    maxlength::Int
    minqual::Int
end

# This converts from PHRED33 to a probability
qual(x::UInt8, offset=UInt8(33)) = 10^-((x - offset) / 10)

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

    return phredmean
end

function passes(record::FASTQ.Record, options::Options)
    if !(options.minlength ≤ length(record.sequence) ≤ options.maxlength)
        return false
    end
    meanqual(record) < options.minqual && return false
    return true
end

function filter(in::IO, out::IO, options::Options)
    reader = FASTQ.Reader(in)
    writer = FASTQ.Writer(out)
    record = FASTQ.Record()
    while !eof(reader)
        read!(reader, record)
        if passes(record, options)
            write(writer, record)
        end
    end
end

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

function maybe_compress(path)
    if endswith(path, ".gz")
        return GzipCompressorStream(open(path, "w"))
    else
        return open(path, "w")
    end
end

@main function main(inpath::String, outpath::String,
              minlength::Int, maxlength::Int,
              minqual::Int)
    options = Options(minlength, maxlength, minqual)
    in = maybe_decompress(inpath)
    out = maybe_compress(outpath)
    filter(in, out, options)
    close(out)
    close(in)
end
