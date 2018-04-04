#!/usr/bin/env julia

# Get model directory and include it in Julia's loading path
#cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/") # Only in ATOM
EDDir = pwd()
push!(LOAD_PATH,EDDir)

srand(123)

# Load Julia & model packages
using Distributions
using JLD
using JuliaDB
using Setworld
using Fileprep
using Organisms

#TODO put them in OrgsRef
const Boltz = 8.62e-5 # Brown & Sibly MTE book chap 2
const aE = 0.65 # Brown & Sibly MTE book chap 2
const plants_gb0 = (10^(10.15))/52 # 10e10.15 is the annual plant biomass production (Ernest et al. 2003) transformed to weekly base
const plants_mb0 = 5.522 #adjustted accordung to 1 death per individual for 1g (MTEpar notebook)
const plants_fb0 = exp(30.0) # fertility rate
const seedmassµ = 0.8
const tK = 273.15 # °C to K converter

"""
    read_initials(simparams, initorgs)
Reads in and stores landscape conditions and organisms from `"landscape_init.in"` and `"organisms.in"` and stores values in composite types.
"""
function read_initials()
    #TODO change simulation input: dictionnary?
    # simparams = Setworld.Simpars()
    # begin
    #     env = open("landscape_init.in")
    #     readline(env); simparams.fxlength = tuple(parse(Int64,readline(env)))
    #     readline(env); simparams.fylength = tuple(parse(Int64,readline(env)))
    #     readline(env); simparams.fmeantemp = tuple(parse(Float64,readline(env)))
    #     readline(env); simparams.ftempsd = tuple(parse(Float64,readline(env)))
    #     readline(env); simparams.fmeanprec = tuple(parse(Float64,readline(env)))
    #     readline(env); simparams.fprecsd = tuple(parse(Float64,readline(env)))
    #     readline(env); simparams.nfrags = parse(Int64,readline(env))
    #     close(env)
    # end
    simparams = Setworld.Simpars()
    simparams.fxlength = (500)
    simparams.fylength = (500)
    simparams.fmeantemp = (20.0)
    simparams.ftempsd = (1.0)
    simparams.fmeanprec = (100.0)
    simparams.fprecsd = 1.0
    simparams.nfrags = 1
    simparams.timesteps = 52
    #verify that: TODO not a real test
    simparams.nfrags == length(simparams.fxlength)

    # initorgs = Organisms.InitOrgs()
    # begin
    #     #TODO check the organisms file format: so far, all fragments get the same sps
    #     orgf = open("organisms.in")
    #     readline(orgf); initorgs.fgroups = tuple(readline(orgf)) #doest need parse for string
    #     readline(orgf); initorgs.sps = tuple(readline(orgf))
    #     readline(orgf); initorgs.init_stage = tuple(readline(orgf))
    #     readline(orgf); initorgs.init_abund = tuple(parse(Int64,readline(orgf)))
    #     #readline(orgf); genotypes = (readline(orgf))
    #     readline(orgf); initorgs.biomassμ = tuple(parse(Float64,readline(orgf)))
    #     readline(orgf); initorgs.biomasssd = tuple(parse(Float64,readline(orgf)))
    #     readline(orgf); initorgs.dispμ = tuple(parse.(Float64,readline(orgf)))
    #     readline(orgf); initorgs.dispshp = tuple(parse.(Float64,readline(orgf)))
    #     readline(orgf); initorgs.radius = tuple(parse(Int64,readline(orgf)))
    #     #check why if parsed into tuple, becomes a float
    #     close(orgf)
    # end
    initorgs = Organisms.InitOrgs()
    initorgs.fgroups = ("wind", "ant")
    initorgs.sps = ("sp1", "sp2")
    initorgs.init_stage = ("a","a")
    initorgs.init_abund = (100,100)
    #genotypes are not initialized with inputs
    initorgs.biomassμ = (100,100)
    initorgs.biomasssd = (1,1)
    initorgs.dispμ = (0,0)
    initorgs.dispshp = (0,0)
    initorgs.radius = (0,0)

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
    outputorgs()
Saves a long format table (`.tsv` file) with the organisms field informations.
"""
function outputorgs(orgs::Array{Organisms.Organism, N} where N, t::Int64)

    #mk dir with simulation parameters identifier

    sep = ","

    output = open(string("EDoutputs/orgsweek",t,".csv"), "w")
    print(output, "id", sep)
    print(output, "location", sep)
    print(output, "species", sep)
    print(output, "stage", sep)
    print(output, "functional_group", sep)
    print(output, "genotype", sep)
    print(output, "dispersal_pars")
    println(output)

    #TODO change writedlm
    for o in 1:length(orgs)
        writedlm(output, [orgs[o].id orgs[o].location orgs[o].sp orgs[o].stage orgs[o].fgroup orgs[o].biomass orgs[o].radius orgs[o].genotype orgs[o].disp])
    end

    close(output)
end

"""
    simulate!()
"""
function simulate()
#   INITIALIZATION
    simparams, initorgs = read_initials()
    mylandscape = landscape_init(simparams)
    orgs = newOrgs(mylandscape, initorgs)

    try
        mkdir("EDoutputs")
    catch
        println("Error in creating output folder (EDoutputs), assuming it already exists")
    end

    simulog = open("EDoutputs/simulog.txt","w")

# MODEL RUN
    for t in 1:simprams.timesteps
        #develop!()
        projvegmass!(mylandscape,orgs,simulog)
        nogrowth = allocate!(mylandscape,orgs,t,aE,Boltz, simulog)
        #TODO check if there is no better way to keep track of individuals that are not growing
        reproduce!(mylandscape,orgs,simulog)
        if (25 <= rem(t, 52) < 37) == 0
            disperse!(mylandscape,orgs,simulog)
            establish!(mylandscape,orgs,simulog)
        end
        survive!(mylandscape, orgs,nogrowth,simulog)
        ## DISTURBANCES
        ## Dynamical landscape change
        # if t #something
        #     function update_landscape!()
        # end
        ## Invasion
        # read_orgs(invasivefile)

        # Output:
        #orgs
        #if rem(t,4) == 0
        outputorgs(orgs,t)
        #end
        #save(string("week",t))
        #network interactions
        #outputnetworks()
        #save(string("week",4*t)) more reasonable interval
    end

    close(simulog)
end

simulate()

#Dictionnary stores functio2nal groups parameters
"""
    fgpars()
Stores functional groups parameters in a dictionnary that is consulted for every function involving organism's simulation.
    !!!! Better than struct in this case because it is possible to write general funcitons that match the strings identifying the fg og the organism and the key in the dictionnary
"""
# function fgparsdict()
#     fgpars = Dict()
#     # parse files in EDDir/functionalgroups and for each file, create an entry in the file with the first 3 letters of the group and the parametr it controls. takes the parameters of a list
#     return fgpars
# end
