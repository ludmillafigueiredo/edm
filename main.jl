#!/usr/bin/env julia

# Get model directory and include it in Julia's loading path
#
#cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/") # Only in ATOM
EDDir = pwd()
push!(LOAD_PATH,EDDir)

srand(123)

# Load Julia & model packages
using ArgParse
using Distributions
#using JLD #saving Julia objects
using JuliaDB #for in/outputs
using DataValues
using Fileprep
using Setworld
using Organisms

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
        required = true
        "--spinput"
        help = "Name of file with species list."
        arg_type = String
        default = abspath(pwd(),"inputs/spinput.csv")
        "--landconfig"
        help = "Name of file with simulation parameters: areas of fragments, mean (and s.d.) temperature, total running time."
        arg_type = String
        default = abspath(pwd(),"inputs/landpars.csv")
        "--timesteps"
        help = "Duration of simulation in weeks."
        arg_type = Int
        default = 52
        "tempts"
        help = "Name of file with weekly temperature  and precipitation time series"
        arg_type = Float
        default = abspath(pwd(),"inputs/envtimeseries.csv")
    end

    return parse_args(sets) # returning a dictionnary of strings is useful because they can passed as keywords to Julia function
end

"""
read_landin(settings)
Reads in and stores landscape conditions and organisms from `"landscape_init.in"` and `"organisms.in"` and stores values in composite types.
"""
function read_landpars(settings::Dict{String,Any})

    #landinputtbl = loadtable(abspath(pwd(),"inputs/landpars.csv"))
    landinputtbl = loadtable(settings["landconfig"])
    temp_tsinput = loadtable(settings["temp_ts"])

    landpars = Setworld.LandPars(Fileprep.areatocell(select(landinputtbl,:areas_m2)),
                                 Fileprep.areatocell(select(landinputtbl,:areas_m2)),
                                 select(temp_tsinput,:meantemp_ts)),
                                 select(temp_tsinput,:sdtemp_ts)),
                                 select(temp_tsinput,:meanprec_ts)),
                                 select(temp_tsinput,:sdprec_ts)),
                                 length(select(landinputtbl, :id))
    return landpars
end

"""
    read_spinput(settings)
Reads in species initial conditions and parameters. Stores tehm in `orgsref`, a structure with parameters names as Dictionnary fields, where species names are the keys to the parameter values.
"""

function read_spinput(settings::Dict{String,Any})

   #spinputtbl = loadtable(abspath(pwd(),"inputs/spinput.csv"))
   spinputtbl = loadtable(settings["spinput"])

   orgsref = OrgsRef(
        Dict(rows(spinputtbl,:sp_id)[i] => rows(spinputtbl,:sp)[i] for i in 1:length(rows(spinputtbl,:sp_id))),
        Array(rows(spinputtbl,:sp_id)),
        Dict(rows(spinputtbl,:sp_id)[i] => rows(spinputtbl,:kernels)[i] for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] => rows(spinputtbl,:biomass_mean)[i] for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] => rows(spinputtbl,:biomass_sd)[i] for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] => rows(spinputtbl,:abund)[i] for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] => rows(spinputtbl,:mean_seed_number)[i] for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] => rows(spinputtbl,:mean_seed_mass)[i] for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] => rows(spinputtbl,:span)[i] for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] => rows(spinputtbl,:max_span)[i] for i in 1:length(rows(spinputtbl,:sp_id)))
    )
    return orgsref
end

"""
    outputorgs(orgs,t,settingsfrgou)
Saves a long format table with the organisms field informations.
"""
function orgstable(orgsref::Organisms.OrgsRef, landpars::Setworld.LandPars, orgs::Array{Organisms.Organism, N} where N, t::Int64, settings::Dict{String,Any})


    if t == 1
        open(string("EDoutputs/",settings["simID"],"/orgsweekly.csv"), "w") do output
            header = append!(["week"], string.(append!(fieldnames(Organism)[1:end-2], ["chrm1" "chrm2" "radius"]))) #for non-diploid Orgs? string.(fieldnames(Organism))
            writedlm(output, reshape(header, 1, length(header)))
        end
    end

    open(string("EDoutputs/",settings["simID"],"/orgsweekly.csv"), "a") do output
        #TODO better extract and arrange the field values
        for o in 1:length(orgs)
            writedlm(output, [t orgs[o].id orgs[o].location orgs[o].sp orgs[o].biomass orgs[o].fgroup orgs[o].stage orgs[o].age orgs[o].reped orgs[o].genotype orgs[o].radius])
        end
    end

    if t == settings["timesteps"]
        open(string("EDoutputs/",settings["simID"],"/simulationID",t), "w") do output
            println(output, "Initial conditions:")
            println(output, dump(orgsref))
            println(output, dump(landpars))
        end
        println("End of simulation")
    end
end

"""
    simulate!()
"""
function simulate()
    #   INITIALIZATION
    settings = parse_commandline()
    #unity test
    println(keys(settings))

    landpars = read_landpars(settings)
    #unity test
    println("Land init stored in object of type $(typeof(landpars))")

    orgsref = read_spinput(settings)
    #unity test
    println("Sp info stored in object of type $(typeof(orgsref))")

    mylandscape = landscape_init(landpars)
    #unity test
    println("Landscape initialized: type $(typeof(mylandscape))")

    orgs = newOrgs(mylandscape, orgsref)
    #unity test
    println("Plants initialized: type $(typeof(orgs))")

    # unity test
    println("Starting simulation")


    try
        mkpath("EDoutputs/$(settings["simID"])")
        println("Output will be written to 'EDoutputs'")
    catch
        println("Overwriting results to existing 'EDoutputs/$(settings["simID"])' folder")
    end

    cd(pwd())

    # OUTPUT SIMULATION SETTINGS
    open(string("EDoutputs/",settings["simID"],"/simID"),"w") do ID
            println(ID, settings)
    end

    # MODEL RUN
    for t in 1:settings["timesteps"]

        println("running week $t")

        # UPDATE TEMPERATURE
        if t != 1
                updatetemp!(mylandscape,t)
        end

        # DISTURBANCE
        #if haskey(settings, "disturb")
        #        # select function according to keyword: loss, frag, temp
        #        if settings["disturb"] == "loss"
        #        elseif settings["disturb"] == "frag"
        #        elseif settings["disturb"] == "temp"
        #        else
        #                error("Please specify one of the disturbance scenarios with `--disturb`:\nhabitat area loss (`loss??),\nhabitat fragmentation (`frag`) or\ntemperature change (`temp`)")
        #                break
        #        end
        #end

        open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                println(sim, "WEEK $t")
        end

        projvegmass!(mylandscape,orgs,settings)

        nogrowth = allocate!(mylandscape,orgs,t,aE,Boltz,settings) #TODO check if there is no better way to keep track of individuals that are not growing

        #juveniles become adults before just before the beggining of spring
        if rem(t, 52) == 11 #juveniles become adults at the beggining of spring (reproductive season)
            develop!(orgs)
        end

        # Plants: adult reproduction and embryos dispersal
        if rem(t, 52) == 24 #reproduction = seed production happens at the end of spring (last week)
            reproduce!(mylandscape,orgs,t, settings, orgsref)
        elseif 25 <= rem(t, 52) < 37  #seed dispersal and germination happen during summer
            disperse!(mylandscape,orgs,settings)
            establish!(mylandscape,orgs,t,settings)
        end

        survive!(mylandscape,orgs,nogrowth,settings) # density-dependent and independent mortality

        # output weekly
        orgstable(orgsref, landpars, orgs,t,settings)

        ## DISTURBANCES
        ##
        ## Dynamical landscape change
        # if t #something
        #     function update_landscape!()
        # end
        ## Invasion
        # read_orgs(invasivefile)

    end

end

simulate()
