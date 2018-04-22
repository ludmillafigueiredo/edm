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
using JuliaDB #for in/outputs
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
        help = "Name of the folder (string type) where outputs will be stored. Default is current time."
        arg_type = String
        default = string(now())
        "--spinfo"
        help = "Name of file with species list."
        arg_type = String
        default = abspath(pwd(),"inputs/spinfo")
        "--landpars"
        help = "Name of file with simulation parameters: areas of fragments, mean (and s.d.) temperature, total running time."
        arg_type = String
        default = abspath(pwd(),"inputs/landpars.csv")
        "--timesteps"
        help = "Duration of simulation (weeks)."
        arg_type = Int64
        default = 52
    end

    return parse_args(sets) # returning a dictionnary of strings is useful because they can passed as keywords to Julia function
end

"""
read_initials(settings)
Reads in and stores landscape conditions and organisms from `"landscape_init.in"` and `"organisms.in"` and stores values in composite types.
"""
function read_initials(settings::Dict{String,Any})

    landin = loadtable(settings["landpars"])
    #TODO if csv cells are not set ot Text, round number are entered as Int64

    landparams = Setworld.Landpars()
    landparams.fxlength = Fileprep.areatocell(select(landin, :areas_m2)) # 5 cm² cells
    landparams.fylength = Fileprep.areatocell(select(landin, :areas_m2))
    landparams.fmeantemp = select(landin, :temp_mean)
    landparams.ftempsd = select(landin, :temp_sd)
    landparams.fmeanprec = select(landin, :precipt_mean)
    landparams.fprecsd = select(landin, :precipt_mean)
    landparams.nfrags == length(landparams.fxlength)
    landparams.timesteps = settings["timesteps"]

    spinfo = loadtable(settings["spinfo"])
    # TODO unique() should be unnecessary once spinfo has a proper format
    initorgs = Organisms.InitOrgs()
    initorgs.fgroups = unique(columns(spinfo, :kernels)) # there is no functional group classification anymore (or yet), so should
    initorgs.sps = unique(columns(spinfo, :sp))
    initorgs.init_stage = "a" #always initializing adults
    initorgs.init_abund = unique(columns(spinfo, :abund))
    #genotypes are not initialized with inputs
    initorgs.biomassμ = unique(columns(spinfo, :biomass_mean))
    initorgs.biomasssd = unique(columns(spinfo, :biomass_sd))
    initorgs.radius = 0

    return landparams, initorgs, spinfo
end

"""
    outputorgs(orgs,t,settingsfrgou)
Saves a long format table with the organisms field informations.
"""
function orgstable(initorgs::Organisms.InitOrgs, landparams::Setworld.Landpars, orgs::Array{Organisms.Organism, N} where N, t::Int64, settings::Dict{String,Any})

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

    if t == seetings["timesteps"]
        open(string("EDoutputs/",settings["simID"],"/simulationID",t,), "w") do output
            println(dump(initorgs))
            println(dump(landparams))
        end
    end
end

"""
    simulate!()
"""
function simulate()
    #   INITIALIZATION
    settings = parse_commandline()
    landparams, initorgs, spinfo = read_initials(settings)
    mylandscape = landscape_init(landparams)
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
    for t in 1:landparams.timesteps

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
