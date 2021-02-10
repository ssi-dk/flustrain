using PackageCompiler 

create_sysimage([:Plots, :FASTX, :BioSequences, :BioAlignments, :CodecZlib, :ErrorTypes, Printf, Serialization],
    sysimage_path="sys_pipe.so", precompile_execution_file="test.jl")
