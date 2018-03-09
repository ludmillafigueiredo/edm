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
#Simulation parameters storage:
mutable struct Simpars #TODO put all landscape.in values in here
    fxlength::Tuple{Int64}
    fylenght::Tuple{Int64}
    fmeantemp::Tuple{Float64}
    ftempsd::Tuple{Float64}
    fmeanprec::Tuple{Float64}
    fprecsd::Tuple{Float64}
    nfrags::Int64
    timesteps::Int64
    Simpars() = new() #necessary
end

# Initial organisms parametrization
mutable struct InitOrgs
    fgroups::Tuple{String}
    sps::Tuple{String}
    init_stage::Tuple{String}
    init_abund::Tuple{String}
    #genotypes #TODO initialie those in functions
    biomassμ::Tuple{Float64}
    biomasssd::Tuple{Float64}
    dispμ::Tuple{Float64}
    dispshp::Tuple{Float64}
    radius::Tuple{Float64}
    InitOrgs() = new() #TODO check if new() is necessary
end

#Dictionnary stores functional groups parameters
"""
    fgpars()
Stores functional groups parameters in a dictionnary that is consulted for every function involving organism's simulation.
    !!!! Better than struct in this case because it is possible to write general funcitons that match the strings identifying the fg og the organism and the key in the dictionnary
"""
function fgparsdict()
    fgpars = Dict()
    # parse files in EDDir/functionalgroups and for each file, create an entry in the file with the first 3 letters of the group and the parametr it controls. takes the parameters of a list
    return fgpars
end

"""
    read_initials(simparams, initorgs)
Reads in and stores landscape conditions and organisms from `"landscape_init.in"` and `"organisms.in"` and stores values in composite types.
"""
function read_initials()
    #TODO dictionnary?
    simparams = Simpars()
    begin
        env = open("landscape_init.in")
        readline(env); simparams.fxlength = parse(Int64,readline(env))
        readline(env); simparams.fylength = parse(Int64,readline(env))
        readline(env); simparams.fmeantemp = parse(Float64,readline(env))
        readline(env); simparams.ftempsd = parse(Float64,readline(env))
        readline(env); simparams.fmeanprec = parse(Float64,readline(env))
        readline(env); simparams.fprecsd = parse(Float64,readline(env))
        readline(env); simparams.n_frags = parse(Int64,readline(env))
        close(env)
    end
    #verify that
    simparams.n_frags == length(simparams.fxlength)

    initorgs = InitOrgs()
    begin
        #TODO check the organisms file format: so far, all fragments get the same sps
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
    return initorgs, simparams
end

# function read_landscape()
#     #read and store initial conditions in dictionnary
#     # for nas linhas do input
#     file = open(readcsv, "landscpinit.csv")
#     map(x,y -> Dict(x => y), file[1,1:end],file[2:end,1:end])
#     landpars = Dict{Any,Float64}(file[1,1:end],file[2:end,1:end])
#
#     #Alternative
#     file = readlines("landscpinit.csv")
#     landpars = Dict()
#     simpars = Simpars()
#     for l in 1:length(file)
#         simpars.???
#         \
#         parse(str, start; greedy=true, raise=true)
#     end
# end

#TODO Store organisms parameters that regulate model run: Might interfere with the initialization values??
mutable struct orgpars
end

"""
"""
function simulate()
#   INITIALIZATION
    read_initials()

    mylandscape= landscape_init(false, simparams)
    orgs_init = newOrgs(mylandscape, initorgs)

    # INITIALIZATION ": multidimensional
    #read_landscape()
    #read_orgs()

# MODEL RUN
    # for t in 1:timesteps
    #     #competition on 2 steps: reflects on compterm
    #     projmass!()
    #     checkcompetition!()
    #     growth!(orgs::Array{Organisms.Organism, N} where N)
    #     reproduction!()
    #     dispersion!()
    #
    #     # Disturbances:
    #     ## Dynamical landscape change
    #     if t #something
    #         function update_landscape!()
    #
    #     ## Invasion
    #     read_orgs(invasivefile)
    #     # Output:
    # end
    return mylandscape, orgs_init
end

simulate()
