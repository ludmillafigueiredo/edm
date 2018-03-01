#!/usr/bin/env julia

# ArgParse propery #TODO

# Get model directory and include in the
cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/")
EDDir = pwd()
push!(LOAD_PATH,EDDir)


global IDcounter = Int64(0)

# Load model packages
using Distributions
using Setworld
#using Organisms

#TODO when reading for more files, consider another input format to handle multiple entries: read as array might be enough

mutable struct simparams #TODO put all landscape.in values in here
    pxlength::Int64
    pylength::Int64
    pmeantemp::Float64
    ptempsd::Float64
    pmeanprec::Float64
    pprecsd::Float64
    n_frags::Int64
    simsettings() = new()
end

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

mutable struct orgparams
    #TODO collect org parameters here
end

# read in which oganisms are going to be simulated
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

mylandscape= landscape_init(false, pxlength, pylength, n_frags, pmeantemp, ptempsd, pmeanprec, pprecsd)
#println("Here is the landscape:")
#write("/mylandscape",mylandscape)

#neighbors = neighborhood(mylandscape)

#TODO use the structs with parameters for it!
orgs_init = newOrgs!(mylandscape,
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

writedlm("mylandscape",dump(mylandscape))

function simulation(
    orgs_init::Array{Organisms.Organism,N} where N,
    neighborhood::Array{Dict, N} where N,
    landscape::Array{Setworld.WorldCell, N} where N
    )
    for t in 1:timesteps
        #competition on 2 steps: reflects on compterm
        projmass!()
        checkcompetition!()
        growth!(orgs::Array{Organisms.Organism, N} where N)
        reproduction!()
        dispersion!()
    end
end
