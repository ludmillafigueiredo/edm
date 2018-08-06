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
using JuliaDB #for in/outputs
using DataValues
using Fileprep
using Setworld
using Organisms

const Boltz = 8.62e-5 # Brown & Sibly MTE book chap 2
const aE = 0.65 # Brown & Sibly MTE book chap 2
#const plants_gb0 = (10^(10.15))/40 # 10e10.15 is the annual plant biomass production (Ernest et al. 2003) transformed to weekly base, with growth not happening during winter (MTEpar notebook)
#const plants_mb0 = 9.902 #adjustted accordung to 1 death per individual for 1g (MTEpar notebook)
#const plants_fb0 = exp(30.0) # fertility rate
#const seedmassµ = 0.8

function parse_commandline()
    sets = ArgParseSettings() #object that will be populated with the arguments by the macro

    @add_arg_table sets begin
        "--simID"
        help = "Name of the folder (string type)
where outputs will be stored."
        arg_type = String
        required = true
        
        "--spinput"
        help = "Name of file with species list."
        arg_type = String
        default = abspath(pwd(),"inputs/species.csv")

        "--tdist"
        help = "Initial trait value distribution"
        arg_type = String
        default = "unif"

        "--insect"
        help = "How to explicitly model insects:
pollination-independent reproduction \"indep\";
equal pollination loss for all species \"equal\"."
        arg_type = String
        default = abspath(pwd(),"inputs/insects.csv")
        
        "--landconfig"
        help = "Name of file with simulation parameters: areas of fragments, mean (and s.d.) temperature, total running time."
        arg_type = String
        default = abspath(pwd(),"inputs/landpars.csv")

        "--disturb"
        help = "Type of environmental disturbance to be implemented: habitat area loss \"loss\", habitat fragmentation \"frag\" or temperature change \"temp\""
        arg_type = String
        default = "none"

        "--timesteps"
        help = "Duration of simulation in weeks."
        arg_type = Int
        default = 52

        "--tout"
        help = "Frequency of output (number of weeks)"
        arg_type = Int
        required = true
        
        "--temp_ts"
        help = "Name of file with weekly temperature  and precipitation time series"
        arg_type = String
        default = abspath(pwd(),"inputs/temp1917_2017.csv")
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
                                 select(temp_tsinput,:meantemp_ts),
                                 select(temp_tsinput,:sdtemp_ts),
                                 select(temp_tsinput,:meanprec_ts),
                                 select(temp_tsinput,:sdprec_ts),
    length(select(landinputtbl, :id)))
    return landpars
end

"""
read_spinput(settings)
Reads in species initial conditions and parameters. Stores tehm in `orgsref`, a structure with parameters names as Dictionnary fields, where species names are the keys to the parameter values.
"""

function read_spinput(settings::Dict{String,Any})

    #spinputtbl = loadtable(abspath(pwd(),"inputs/species.csv"))
    spinputtbl = loadtable(settings["spinput"])

    orgsref = OrgsRef(Array(rows(spinputtbl,:sp_id)),
                      Dict(rows(spinputtbl,:sp_id)[i] =>
                           rows(spinputtbl,:kernel)[i]
                           for i in 1:length(rows(spinputtbl,:sp_id))),
                      Dict(rows(spinputtbl,:sp_id)[i] =>
                           rows(spinputtbl,:kernel_sd)[i]
                           for i in 1:length(rows(spinputtbl,:sp_id))),
                      Dict(rows(spinputtbl,:sp_id)[i] =>
                           rows(spinputtbl,:e_mu)[i]
                           for i in 1:length(rows(spinputtbl,:sp_id))),
                      Dict(rows(spinputtbl,:sp_id)[i] =>
                           rows(spinputtbl,:e_sd)[i]
                           for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0g)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0g_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0em)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0em_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0am)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0am_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0jg)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0jg_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0ag)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0ag_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:sestra)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:dyad)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:floron)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:floron_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:floroff)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:floroff_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:seedon)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:seedon_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:seedoff)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:seedoff_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:max_mass)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:min_mass)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:min_mass_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:max_span)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:max_span_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:mass_mu)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:mass_sd)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:abund)[i]
         for i in 1:length(rows(spinputtbl,:sp_id)))
    )
    return orgsref
end

"""
implicit_insect(settings)
Reads how insects are going to be implicitly simulated.
"""
function implicit_insect(settings::Dict{String,Any})

    insectsinput = loadtable(settings["insect"])
    
    interaction = select(insectsinput, :interaction)[1]
    scen = select(insectsinput, :scen)[1] # pollination scenario
    tp = select(insectsinput, :tp)[1] # time of perturbation
    remain = select(insectsinput, :remain)[1] # remain proportion of insects right after pertubation
    regime = select(insectsinput, :regime)[1] # regime of loss
    
    return interaction, scen, tp, remain, regime
end

"""
outputorgs(orgs,t,settingsfrgou)
Saves a long format table with the organisms field informations.
"""
function orgstable(orgsref::Organisms.OrgsRef, landpars::Setworld.LandPars, orgs::Array{Organisms.Organism, N} where N,landscape::Array{Dict{Any,Any}}, t::Int64, settings::Dict{String,Any})

    if t == 1
        header = hcat(["week"],
                      reshape(string.(fieldnames(Organism)[1:3]),1,3),
                      ["veg" "repr"],
                      reshape(string.(fieldnames(Organism)[5:end-2]),1,length(fieldnames(Organism)[5:end-2])),
                      ["chrm1" "chrm2" "radius"])#for non-diploid Orgs? string.(fieldnames(Organism))
        open(string("EDoutputs/",settings["simID"],"/orgsweekly.csv"), "w") do output
            writedlm(output, header) #reshape(header, 1, length(header)))
        end
        open(string("EDoutputs/", settings["simID"], "/landscape.csv"), "w") do landoutput
            println(landoutput, landscape[1,1:end,1])
        end
    end

    if rem(t,settings["tout"]) == 0 # output bi-monthly
        for o in 1:length(orgs)
            open(string("EDoutputs/",settings["simID"],"/orgsweekly.csv"), "a") do output
                writedlm(output, hcat(t,
                                      orgs[o].id,
                                      orgs[o].location,
                                      orgs[o].sp,
                                      orgs[o].mass["veg"],
                                      orgs[o].mass["repr"],
                                      orgs[o].stage,
                                      orgs[o].age,
                                      orgs[o].mated,
                orgs[o].genotype,
                orgs[o].radius))
            end
        end
        #TODO output orgsref and ladpars
        open(string("EDoutputs/", settings["simID"], "/landscape.csv"), "a") do landoutput
            println(landoutput,landscape[1,1:end,1])
        end
    end
end
"""
disturb!()

"""
function disturb!(landscape::Array{Dict{Any,Any}}, landavail::Array{Bool,N} where N, orgs::Array{Organisms.Organism, N} where N, t::Int64, settings::Dict{String,Any})
    # read in the disturbance file, if not done so yet
    #if !isdefined(:disturbtbl)
        # select file according to keyword: loss, frag, temp
        if settings["disturb"] == "loss"
            disturbtbl = loadtable(abspath(pwd(),"inputs/arealoss.csv"))
            tdist = select(disturbtbl,:time)
        elseif settings["disturb"] == "frag"
            disturbtbl = loadtable(abspath(pwd(),"inputs/fragmentation.csv"))
            tdist = select(disturbtbl,:time)
        elseif settings["disturb"] == "temp"
            #continue
            println("Temperature change is simulated with
                    the provided temperature file.")
        else
            error("Please specify one of the disturbance scenarios with `--disturb`:
                  \n\"none\" if no disturbance should be simulated,
                  \n\"loss\" for habitat area loss,
                  \n\"frag\" for habitat fragmentation,
                  \n\"temp\" for temperature change.")
        end
    #else
        if t in tdist
            if settings["disturb"] == "loss"
                loss = select(filter(x -> x.time == t,disturbtbl),
                              :proportion)[1] # select returns an array
                Setworld.destroyarea!(landavail,loss,settings)
                Organisms.destroyorgs!(orgs,landavail,settings)
            elseif settings["disturb"] == "frag"
                # fragment!(mylandscape,orgs)
                #while fragment is not implemented
            end
        end
    #end
    #return disturbtbl, tdist
end

"""
simulate!()
"""
function simulate()
    # INITIALIZATION
    # Read in command line arguments
    settings = parse_commandline()
    # unity test
    println(keys(settings))

    # Store landscape configuration
    landpars = read_landpars(settings)
    # unity test
    println("Land init stored in object of type $(typeof(landpars))")

    # Store species information
    orgsref = read_spinput(settings)
    # unity test
    println("Sp info stored in object of type $(typeof(orgsref))")

    # Set insects implicit simulation
    interaction, scen, tp, remain, regime = implicit_insect(settings)
    
    # Create landscape
    mylandscape, landavail = landscape_init(landpars)
    # unity test
    println("Landscape initialized: type $(typeof(mylandscape))")

    # Initialize individual tagger: It tags all individuals ever created. It does not change if individuals die, so there is no risk of re-use of a tag.
    id_counter = 0
    
    # Create initial individuals
    orgs, id_counter = newOrgs!(landavail, orgsref, id_counter, settings["tdist"])
    
    println("Plants initialized: type $(typeof(orgs))")

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
    # START ID SIMULATION LOG FILE
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"w") do sim
        println(sim,string("Simulation: ",settings["simID"],now()))
    end

    # MODEL RUN
    for t in 1:settings["timesteps"]

        println("running week $t")
        # unity test
        open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            println(sim, "WEEK $t")
        end

        #UPDATE LANDSCAPE: weekly temperature and projected biomass
        T = updateenv!(mylandscape, t, landpars)
        
        # DISTURBANCE
        if settings["disturb"] != "none"  
           disturb!(mylandscape,landavail,orgs,t,settings)
        end

        # LIFE CYCLES
        projvegmass!(mylandscape,orgs, settings)

        nogrowth = allocate!(mylandscape,orgs,t,aE,Boltz,settings,orgsref, T)

        develop!(orgs,orgsref)

        mate!(orgs,t,settings,scen,tp,remain,regime, 1)

        id_counter = mkoffspring!(orgs,t,settings,orgsref,id_counter)

        seedsi = release!(orgs,t,settings,orgsref)

        disperse!(landavail,seedsi,orgs,t,settings,orgsref)

        establish!(mylandscape,orgs,t,settings,orgsref)
        
        survive!(orgs,nogrowth,t,settings,orgsref,T)

        shedd!(orgs,orgsref,t)

        # OUTPUT
        orgstable(orgsref,landpars,orgs,mylandscape,t,settings)
    end
end
simulate()
