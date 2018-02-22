#!/usr/bin/env julia

# ArgParse propery #TODO

# Get model directory and include in the
cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/")
EDDir = pwd()
push!(LOAD_PATH,EDDir)

# Load modules:
# General
global ED_SEED =srand(123456789) # set seed for randomisation

# Load model packages
using Distributions
using Setworld
using Organisms

#TODO when reading for more files, consider another input format to handle multiple entries: read as array might be enough
begin
    env = open("landscape_init.in")
    readline(env); pxlength = parse(Int64,readline(env))
    readline(env); pylength = parse(Int64,readline(env))
    readline(env); pmeantemp = parse(Float64,readline(env))
    readline(env); ptempsd = parse(Float64,readline(env))
    readline(env); pmeanprec = parse(Float64,readline(env))
    readline(env); pprecsd = parse(Float64,readline(env))
    readline(env); n_frags = parse(Int64,readline(env))
    close(env)
end

# read in which oganisms are going to be simulated
# TODO stores it in some kind of
begin
    orgf = open("organisms.in")
    readline(orgf); fgroups = [readline(orgf)] #doest need parse for string
    readline(orgf); sps = [readline(orgf)]
    readline(orgf); init_stage = [readline(orgf)]
    readline(orgf); init_abund = [parse(Int64,readline(orgf))]
    readline(orgf); genotypes = [readline(orgf)]
    readline(orgf); biomassμ = [parse(Float64,readline(orgf))]
    readline(orgf); biomasssd = [parse(Float64,readline(orgf))]
    readline(orgf); dispμ = [parse.(Float64,readline(orgf))]
    readline(orgf); dispsd = [parse.(Float64,readline(orgf))]
    readline(orgf); radius = [parse(Int64,readline(orgf))]
    close(orgf)
end

# Initialize individuals #TODO do proper simulation initialization
# TODO worth keeping this function? just 3 lines

mylandscape = landscape_init(false, pxlength, pylength, n_frags, pmeantemp, ptempsd, pmeanprec, pprecsd)
println("Here is the landscape:")
write("/mylandscape",mylandscape)


newOrgs!(mylandscape,
fgroups,
sps,
init_stage,
init_abund,
biomassμ,
biomasssd,
genotypes,
dispμ,
dispsd,
radius,
IDcounter)
