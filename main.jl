#!/usr/bin/env julia

# Get model directory and include it
cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/") # Only in ATOM
EDDir = pwd()
push!(LOAD_PATH,EDDir)

global IDcounter = Int64(0)

# Load Julia & model packages
using Distributions
using JLD
using Setworld
#using Organisms

#Simulation parameters storage:
mutable struct Simpars #TODO put all landscape.in values in here
    fxlength::Tuple{Int64}
    fylength::Tuple{Int64}
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
    init_abund::Tuple{Int64}
    #genotypes #TODO initialie those in functions
    biomassμ::Tuple{Float64}
    biomasssd::Tuple{Float64}
    dispμ::Tuple{Float64}
    dispshp::Tuple{Float64}
    radius::Tuple{Int64}
    InitOrgs() = new() #TODO check if new() is necessary
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
        readline(env); simparams.fxlength = tuple(parse(Int64,readline(env)))
        readline(env); simparams.fylength = tuple(parse(Int64,readline(env)))
        readline(env); simparams.fmeantemp = tuple(parse(Float64,readline(env)))
        readline(env); simparams.ftempsd = tuple(parse(Float64,readline(env)))
        readline(env); simparams.fmeanprec = tuple(parse(Float64,readline(env)))
        readline(env); simparams.fprecsd = tuple(parse(Float64,readline(env)))
        readline(env); simparams.nfrags = parse(Int64,readline(env))
        close(env)
    end
    #verify that
    simparams.nfrags == length(simparams.fxlength)

    initorgs = InitOrgs()
    begin
        #TODO check the organisms file format: so far, all fragments get the same sps
        orgf = open("organisms.in")
        readline(orgf); initorgs.fgroups = tuple(readline(orgf)) #doest need parse for string
        readline(orgf); initorgs.sps = tuple(readline(orgf))
        readline(orgf); initorgs.init_stage = tuple(readline(orgf))
        readline(orgf); initorgs.init_abund = tuple(parse(Int64,readline(orgf)))
        #readline(orgf); genotypes = (readline(orgf))
        readline(orgf); initorgs.biomassμ = tuple(parse(Float64,readline(orgf)))
        readline(orgf); initorgs.biomasssd = tuple(parse(Float64,readline(orgf)))
        readline(orgf); initorgs.dispμ = tuple(parse.(Float64,readline(orgf)))
        readline(orgf); initorgs.dispshp = tuple(parse.(Float64,readline(orgf)))
        readline(orgf); initorgs.radius = tuple(parse(Int64,readline(orgf)))
        #check why if parsed into tuple, becomes a float
        close(orgf)
    end
    return simparams, initorgs
end

# function read_landscape()
#     #read and store initial conditions in dictionnary
#     # for nas linhas do input
#     file = open(readcsv, "landscpinit.csv")
#     map(x,y -> Dict(x => y), file(1,1:end),file(2:end,1:end))
#     landpars = Dict{Any,Float64}(file(1,1:end),file(2:end,1:end))
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


"""
"""
function simulate(timesteps)
#   INITIALIZATION
    simparams, initorgs = read_initials()

    mylandscape = landscape_init(simparams)
    orgs = newOrgs(mylandscape, initorgs)

    # INITIALIZATION ": multidimensional
    #read_landscape()
    #read_orgs()

# MODEL RUN
    for t in 1:timesteps
        #competition on 2 steps: reflects on compterm
        projvegmass!(mylandscape,orgs)
        allocate!()
        reproduce!()
        disperse!()
        survive!()
        # Disturbances:
        ## Dynamical landscape change
        # if t #something
        #     function update_landscape!()
        # end
        ## Invasion
        # read_orgs(invasivefile)

        # Output:
    end
    return mylandscape, orgs_init
end

simulate()

#TODO Store organisms parameters that regulate model run: Might interfere with the initialization values??
mutable struct orgpars
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
