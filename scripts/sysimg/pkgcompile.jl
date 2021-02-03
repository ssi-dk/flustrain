using PackageCompiler 

create_sysimage([:Plots, :FASTX, :BioSequences, :SHA, :CodecZlib, :Comonicon], sysimage_path="sys_pipe.so", precompile_execution_file="test.jl")
