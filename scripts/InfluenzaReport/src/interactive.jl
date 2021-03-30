BLASTN = "/Users/jakobnissen/miniconda3/bin/blastn"
REFDIR = "/Users/jakobnissen/Documents/ssi/projects/flupipe/ref/seqs"

"""
    manual_check(::Segment, ::LongDNASeq; refdir=[preset])

Checks a segment sequence by BLASTing it to a collection of references.
Returns a (ReferenceAssembly, Vector{String}) result, with the strings being
the lines would be in a flupipe report.

Note that this approach is fairly inefficient, and should probably not be used in
a loop unless you want your CPU to keep you warm in the winter.
"""
function manual_check(segment::Segment, seqs, refdir::String=REFDIR)
    seqlist = collect(seqs)

    # Write to file
    filename = tempname()
    open(filename, "w") do file
        for (i, seq) in enumerate(seqlist)
            println(file, ">$i\n", seq)
        end
    end

    # BLAST to ref
    outfile = tempname()
    subject = joinpath(refdir, "$segment.fna")
    command = `$BLASTN -query $filename -subject $subject -outfmt "6 qacc sacc bitscore qlen length pident"`
    run(pipeline(command, stdout=outfile))

    best_hits = open(outfile) do file
        x = parse_blastout(file, 0.9)
        len = length(x)
        d = Dict(parse(Int, k) => expect(v, "No good hits found") for (k,v) in x)
        @assert length(d) == len
        return d
    end

    # Extract best hits, get that hit from jls file
    references = load_references(segment, refdir, Set(values(best_hits)))
    
    # Create ReferenceAssemblies
    result = Vector{ReferenceAssembly}()
    for (index, seq) in enumerate(seqlist)
        hit_name = best_hits[index]
        assembly = Assembly(segment, falses(length(seq)), seq, hit_name)
        push!(result, ReferenceAssembly(references[hit_name], assembly))
    end

    return result
end

function parse_blastout(io::IO, lenratio::Real)::Dict{String, Option{String}}
    fields = Vector{SubString{String}}(undef, 6)
    hits = eachline(io) |> Map(strip) |> Filter(!isempty) |> Map() do line
        split!(fields, line, UInt8('\t'))
        bitscore = parse(Float64, fields[3])
        qlen = parse(UInt, fields[4])
        len = parse(UInt, fields[5])
        ident = parse(Float64, fields[6]) / 100
        (String(first(fields)), String(fields[2]), bitscore, qlen, len, ident)
    end |> collect
    
    # Group by query sequence
    byquery = Dict{String, Vector}()
    for hit in hits
        push!(get!(Vector{Any}, byquery, first(hit)), hit[2:end])
    end

    foreach(values(byquery)) do hits
        sort!(hits, by=x -> x[2], rev=true)
    end

    # Get best hits
    result = Dict{String, Option{String}}()
    for (query, hits) in byquery
        for (subject, bitscore, qlen, len, ident) in hits
            if (len / qlen) ≥ lenratio && ident ≥ 0.9
                result[query] = some(subject)
                break
            end
        end
        haskey(result, query) || (result[query] = none(String))
    end
    return result
end
