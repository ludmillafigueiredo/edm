#!/usr/bin/env julia

# ArgParse propery #TODO

# Get model directory and include in the
cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/")
EDDir = pwd()
push!(LOAD_PATH,EDDir)

# Load modules:
# General
global GLOBAL_SEED =srand(123) # set seed for randomisation

# World building
using setworld

#TODO when reading for more files, considering transforming it in function to handle multiple entries
begin
    env = open("temp_prec.in")
    readline(env); pxlength = parse(Int64,readline(env))
    readline(env); pylength = parse(Int64,readline(env))
    readline(env); pmeantemp = parse(Float64,readline(env))
    readline(env); ptempsd = parse(Float64,readline(env))
    readline(env); pmeanprec = parse(Float64,readline(env))
    readline(env); pprecsd = parse(Float64,readline(env))
    readline(env); n_frags = parse(Int64,readline(env))
    close(env)
end
println("Environmental input read.")

mylandscape = landscape_init(false, pxlength, pylength, n_frags, pmeantemp, ptempsd, pmeanprec, pprecsd)

println("Here is the landscape:")
#TODO write("/mylandscape",mylandscape)
