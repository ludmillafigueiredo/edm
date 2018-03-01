#!/usr/bin/env julia

# Get model directory and include it
cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/") # Only in ATOM
EDDir = pwd()
push!(LOAD_PATH,EDDir)

global IDcounter = Int64(0)

# Load Julia & model packages
using Distributions
using Setworld
#using Organisms

#TODO julia struct style: caps
#Parameter storage:
mutable struct simpars #TODO put all landscape.in values in here
    pxlength::Int64
    pylength::Int64
    pmeantemp::Float64
    ptempsd::Float64
    pmeanprec::Float64
    pprecsd::Float64
    n_frags::Int64
    simpars() = new()
end

mutable struct initorgvalues
    fgroups
    sps
    init_stage
    init_abund
    #genotypes #TODO initialie those in functions
    biomassμ
    biomasssd
    dispμ
    dispshp
    radius
    initorgvalues() = new() #TODO check if new() is necessary
end

"""
    read_initials(simparams, initorgs)
Reads in and stores landscape conditions and organisms from `"landscape_init.in"` and `"organisms.in"` and stores values in composite types.
"""
function read_initials!(simparams, initorgs)
    #TODO dictionnary
begin
    env = open("landscape_init.in")
    readline(env); simparams.pxlength = parse(Int64,readline(env))
    readline(env); simparams.pylength = parse(Int64,readline(env))
    readline(env); simparams.pmeantemp = parse(Float64,readline(env))
    readline(env); simparams.ptempsd = parse(Float64,readline(env))
    readline(env); simparams.pmeanprec = parse(Float64,readline(env))
    readline(env); simparams.pprecsd = parse(Float64,readline(env))
    readline(env); simparams.n_frags = parse(Int64,readline(env))
    close(env)
end

begin
    orgf = open("organisms.in")
    readline(orgf); initorgs.fgroups = [readline(orgf)] #doest need parse for string
    readline(orgf); initorgs.sps = [readline(orgf)]
    readline(orgf); initorgs.init_stage = [readline(orgf)]
    readline(orgf); initorgs.init_abund = [parse(Int64,readline(orgf))]
    #readline(orgf); genotypes = [readline(orgf)]
    readline(orgf); initorgs.biomassμ = [parse(Float64,readline(orgf))]
    readline(orgf); initorgs.biomasssd = [parse(Float64,readline(orgf))]
    readline(orgf); initorgs.dispμ = [parse.(Float64,readline(orgf))]
    readline(orgf); initorgs.dispshp = [parse.(Float64,readline(orgf))]
    readline(orgf); initorgs.radius = [parse(Int64,readline(orgf))]
    close(orgf)
end
end


function read_landscape()
    try
        #read and store initial conditions in dictionnary
    catch
        #read new conditions and change dictionnary
        # there should be some control for when this change happens
    end
end

#TODO Store organisms parameters that regulate model run: Might interfere with the initialization values??
mutable struct orgpars
end

"""
"""
function simulate()
#   INITIALIZATION
    simparams = simpars()
    initorgs = initorgvalues()
    read_initials!(simparams,initorgs)

    mylandscape= landscape_init(false, simparams)
    orgs_init = newOrgs(mylandscape, initorgs)

# MODEL RUN
    for t in 1:timesteps
        #competition on 2 steps: reflects on compterm
        projmass!()
        checkcompetition!()
        growth!(orgs::Array{Organisms.Organism, N} where N)
        reproduction!()
        dispersion!()
    end
end

writedlm("mylandscape",dump(mylandscape))
