#!/usr/bin/env julia

# Get model directory and include it in Julia's loading path
#cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/") # Only in ATOM
workspace()
EDDir = pwd()
push!(LOAD_PATH,EDDir)

srand(123)

# Load Julia & model packages
using ArgParse
using Distributions
#using JLD #saving Julia objects
using JuliaDB #for outputs
using Setworld
using Fileprep
using Organisms

#TODO put them in OrgsRef
const Boltz = 8.62e-5 # Brown & Sibly MTE book chap 2
const aE = 0.65 # Brown & Sibly MTE book chap 2
const plants_gb0 = (10^(10.15))/40 # 10e10.15 is the annual plant biomass production (Ernest et al. 2003) transformed to weekly base, with growth not happening during winter (MTEpar notebook)
const plants_mb0 = 1.5029220413821088e11 #adjustted accordung to 1 death per individual for 1g (MTEpar notebook)
const plants_fb0 = exp(30.0) # fertility rate
const seedmassµ = 0.8

function parse_commandline()
    sets = ArgParseSettings() #object that will be populated with the arguments by the macro

    @add_arg_table sets begin
        "--simID"
        help = "Name of the folder (string type) where outputs are stored. Default is current time."
        arg_type = String
        default = string(now())
        "--spfile"
        help = "Name of file with species list."
        arg_type = String
        default = string(pwd(),"splist.txt")
    end

    return parse_args(sets) # returning a dictionnary of strings is useful because they can passed as keywords to Julia function
end

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
    simparams.fxlength = tuple(50) # 5 cm² cells
    simparams.fylength = tuple(50)
    simparams.fmeantemp = tuple(20.0)
    simparams.ftempsd = tuple(1.0)
    simparams.fmeanprec = tuple(100.0)
    simparams.fprecsd = tuple(1.0)
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
    initorgs.fgroups = tuple("wind", "ant")
    initorgs.sps = tuple("sp1", "sp2")
    initorgs.init_stage = tuple("a","a")
    initorgs.init_abund = tuple(1,1)
    #genotypes are not initialized with inputs
    initorgs.biomassμ = tuple(100,100)
    initorgs.biomasssd = tuple(1,1)
    initorgs.dispμ = tuple(0,0)
    initorgs.dispshp = tuple(0,0)
    initorgs.radius = tuple(0,0)

    return simparams, initorgs
end

# """
#     input(orgsfilepath,landfilepath)
# Reads in files describing initial species pool, along with functional group, life stage, abundance and weight(mean and ⨦sd).
# """
#
# function input()
#
#     species = loadtable("species.csv")
#     #fg1pars = Dict("sps" => ["sp1", "sp2", "sp3"], "mean" => 100, "sd" => 1)
#     #fg1pars = Dict("sps" => ["sp1", "sp2", "sp3"], "mean" => 100, "sd" => 1)
# end

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
    outputorgs(orgs,t,settingsfrgou)
Saves a long format table with the organisms field informations.
"""
function orgstable(orgs::Array{Organisms.Organism, N} where N, t::Int64, settings::Dict{String,Any})

    #mk dir with simulation parameters identifier

    sep = ','

    if t == 1
        open(string("EDoutputs/",settings["simID"],"/orgsweek",t,".csv"), "a") do output
            header = append!(["week"], string.(fieldnames(Organisms)))
            writedlm(output, reshape(header, 1, length(header)), sep)
        end
    end

    open(string("EDoutputs/",settings["simID"],"/orgsweek",t,".csv"), "a") do output
        #TODO better extract and arrange the variables
        for o in 1:length(orgs)
            writedlm(output, [t orgs[o].id orgs[o].sp orgs[o].stage orgs[o].fgroup orgs[o].location sum(values(orgs[o].biomass)) orgs[o].radius orgs[o].genotype orgs[o].disp], sep)
        end
    end

end

"""
    simulate!()
"""
function simulate()
    #   INITIALIZATION
    settings = parse_commandline()
    simparams, initorgs = read_initials()
    mylandscape = landscape_init(simparams)
    orgs = newOrgs(mylandscape, initorgs)

    # unity test
    println("Starting simulation")

    try
        mkpath("EDoutputs/$(settings["simID"])")
        println("Output will be written to 'EDoutputs'")
    catch
        println("Writing results to existing 'EDoutputs/$(settings["simID"])' folder")
    end

    cd(pwd())

    # MODEL RUN
    for t in 1:simparams.timesteps

        println("running week $t")

        #unity testing
        open(string("EDoutputs/",settings["simID"],"/orgsweek",t,".csv"), "a") do sim
            println(sim,"WEEK ",t)
        end

        projvegmass!(mylandscape,orgs,settings)

        nogrowth = allocate!(mylandscape,orgs,t,aE,Boltz,settings) #TODO check if there is no better way to keep track of individuals that are not growing

        #juveniles become adults before just before the beggining of spring
        if rem(t, 52) == 11 #juveniles become adults at the beggining of spring (reproductive season)
            develop!(orgs)
        end

        # Plants: adult reproduction and embryos dispersal
        if rem(t, 52) == 24 #reproduction = seed production happens at the end of spring (last week)
            reproduce!(mylandscape,orgs,t, settings)
        elseif 25 <= rem(t, 52) < 37  #seed dispersal and germination happen during summer
            disperse!(mylandscape,orgs,settings)
            establish!(mylandscape,orgs,t,settings)
        end

        survive!(mylandscape,orgs,nogrowth,settings) # density-dependent and independent mortality

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
        orgstable(orgs,t,settings)
        #end
        #save(string("week",t))
        #network interactions
        #outputnetworks()
        #save(string("week",4*t)) more reasonable interval
    end
end

simulate()

#Dictionnary stores functio2nal groups parameters
# """
# fgpars()
# Stores functional groups parameters in a dictionnary that is consulted for every function involving organism's simulation.
# !!!! Better than struct in this case because it is possible to write general funcitons that match the strings identifying the fg og the organism and the key in the dictionnary
# """
    # function fgparsdict()
    #     fgpars = Dict()
    #     # parse files in EDDir/functionalgroups and for each file, create an entry in the file with the first 3 letters of the group and the parametr it controls. takes the parameters of a list
    #     return fgpars
    # end
