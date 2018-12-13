#!/usr/bin/env julia

# Get model directory and include it in Julia's loading path
#
#cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/") # Only in ATOM
EDDir = pwd()
push!(LOAD_PATH,EDDir)

# Load Julia & model packages
using ArgParse
using Distributions
using JuliaDB #for in/outputs
using DataValues
using Fileprep
using Setworld
using Organisms

#const Boltz = 1.38064852e-23 # Alternatively:
const Boltz = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
#const aE = 1e-19
const aE = 0.63
global K = 0.0
global cK = 0.0
global Ah = 0.0
global AT = 0.0
global nogrowth = Int64[]

function parse_commandline()
    sets = ArgParseSettings() #object that will be populated with the arguments by the macro

    @add_arg_table sets begin
        "--simID"
        help = "Name of the folder where outputs will be stored."
        arg_type = String
        required = true
        
	"--rseed"
	help = "Seed for RNG"
	arg_type = Int
	required = true

        "--outputat"
        help = "Name of directory where output should be written ."
        arg_type = String
        default = abspath(pwd(),"EDoutputs")
	
	"--spinput"
        help = "Name of file with species list."
        arg_type = String
        default = abspath(pwd(),"inputs/species.csv")

        "--competition"
        help = "Type of competition: \"individual\", based on FON or \"capacity\", based on the landscape carrying capacity"
        arg_type = String
        default = "capacity"

        "--insect"
        help = "How to explicitly model insects:
            pollination-independent reproduction \"indep\";
            equal pollination loss for all species \"equal\"."
        arg_type = String
        default = abspath(pwd(),"inputs/insects.csv")
        
        "--landconfig"
        help = "Name of file with landscape size values: areas of fragments, mean (and s.d.) temperature."
        arg_type = String
        default = abspath(pwd(),"inputs/landpars.csv")

        "--connects"
        help = "Connectivity matrix of distances between fragments."
        arg_type = String
        default = abspath(pwd(),"inputs/connects")
        
        "--disturb"
        help = "Type of environmental disturbance to be implemented: habitat area loss \"loss\", habitat fragmentation \"frag\" or temperature change \"temp\""
        arg_type = String
        default = "none"
        
        "--areafile"
        arg_type = String
        default = abspath(pwd(),"inputs/arealoss.csv")

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

        "--timemsg"
        help = "Output timing to terminal, as well as simulog."
        arg_type = Bool
        default = false
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
         rows(spinputtbl,:e_long)[i]
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
         rows(spinputtbl,:b0jm)[i]
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:b0jm_sd)[i]
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
         rows(spinputtbl,:sestra)[i] == "true"
         for i in 1:length(rows(spinputtbl,:sp_id))),
    Dict(rows(spinputtbl,:sp_id)[i] =>
         rows(spinputtbl,:max_seedn)[i]
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
    regime = select(insectsinput, :regime)[1] # regime of loss
    
    return interaction, scen, tp, regime
end

"""
outputorgs(orgs,t,settingsfrgou)
Saves a long format table with the organisms field informations.
"""
function orgstable(orgsref::Organisms.OrgsRef, landpars::Setworld.LandPars, orgs::Array{Organisms.Organism,1},landscape::Array{Dict{String,Float64},2}, t::Int64, settings::Dict{String,Any})

    outorgs = find(x -> x.stage != "e", orgs)
    
    if t == 1
        header = hcat(["week"],
                      reshape(string.(fieldnames(Organism)[1:4]),1,4),
                      ["veg" "repr"],
                      reshape(string.(fieldnames(Organism)[6:end]),1,length(fieldnames(Organism)[6:end])))#for non-diploid Orgs? string.(fieldnames(Organism))
        open(abspath(joinpath(settings["outputat"],settings["simID"],"orgsweekly.txt")), "w") do output
            writedlm(output, header) #reshape(header, 1, length(header)))
        end
        
        for o in outorgs
            open(abspath(joinpath(settings["outputat"],settings["simID"],"orgsweekly.txt")), "a") do output
                writedlm(output, hcat(t,
                                      orgs[o].id,
                                      orgs[o].stage,
                                      orgs[o].location,
                                      orgs[o].sp,
                                      orgs[o].mass["veg"],
                                      orgs[o].mass["repr"],
                                      orgs[o].kernel,
                                      orgs[o].e_mu,
                orgs[o].seedbank,
                orgs[o].b0g,
                orgs[o].b0em,
                orgs[o].b0jm,
                orgs[o].b0am,
                orgs[o].b0jg,
                orgs[o].b0ag,
                orgs[o].sestra,
                orgs[o].floron,
                orgs[o].floroff,
                orgs[o].wseedn,
                orgs[o].seedon,
                orgs[o].seedoff,
                orgs[o].max_mass,
                orgs[o].max_span,
                orgs[o].age,
                orgs[o].mated))
            end
        end
        
        open(abspath(joinpath(settings["outputat"], settings["simID"], "landscape.csv")), "w") do landoutput
            println(landoutput, landscape[1,1:end,1])
        end
    end

    if rem(t,settings["tout"]) == 0 
        for o in outorgs
            open(abspath(joinpath(settings["outputat"],settings["simID"],"orgsweekly.txt")), "a") do output
                writedlm(output, hcat(t,
                                      orgs[o].id,
                                      orgs[o].stage,
                                      orgs[o].location,
                                      orgs[o].sp,
                                      orgs[o].mass["veg"],
                                      orgs[o].mass["repr"],
                                      orgs[o].kernel,
                                      orgs[o].e_mu,
                orgs[o].seedbank,
                orgs[o].b0g,
                orgs[o].b0em,
                orgs[o].b0jm,
                orgs[o].b0am,
                orgs[o].b0jg,
                orgs[o].b0ag,
                orgs[o].sestra,
                orgs[o].floron,
                orgs[o].floroff,
                orgs[o].wseedn,
                orgs[o].seedon,
                orgs[o].seedoff,
                orgs[o].max_mass,
                orgs[o].max_span,
                orgs[o].age,
                orgs[o].mated))
            end
        end
        #TODO output orgsref and ladpars
        open(abspath(joinpath(settings["outputat"], settings["simID"], "landscape.csv")), "a") do landoutput
            println(landoutput,landscape[1,1:end,1])
        end
    end
end
"""
disturb!()

"""
function disturb!(landscape::Array{Dict{String,Float64},2}, landavail::Array{Bool,2}, orgs::Array{Organisms.Organism,1}, t::Int64, settings::Dict{String,Any})
    # read in the disturbance file, if not done so yet
    #if !isdefined(:disturbtbl)
    # select file according to keyword: loss, frag, temp
    if settings["disturb"] != "none"
        if settings["disturb"] == "loss"
            disturbtbl = loadtable(settings["areafile"])
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
        return tdist
    end
end

function losschange(landavail::Array{Bool,2}, settings::Dict{String,Any}, t::Int64, tdist::Any)
    
    if t == 1 || (settings["disturb"] == "loss" && t in [(tdist-1) tdist (tdist-1)])
        if t == 1
            open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"w") do sim
                println(sim, "week\tK\tcK")
            end          
        end
        global K = (1/100)*(length(find(x -> x == true, landavail))*25) #x tons/ha = x.100g/1m²
        global cK = K/length(find(x -> x == true, landavail))
        
        open(abspath(joinpath(settings["outputat"],settings["simID"],"landlog.txt")),"a") do sim
            writedlm(sim,[t K cK])
        end
    end
    
end

function fragchange(landavail::Array{Bool,2}, settings::Dict{String,Any}, t::Int64, connects::Array{Float64,2}, tdist::Any)
    if t == 1 || (settings["disturb"] == "frag" && t == tdist) 
        global Ah = length(find(x -> x == true, landavail))*25*0.0001
        global AT = sort(reshape(connects,prod(size(connects))), rev = true) |> (y -> length(y) > 1 ? prod(y[1:2]) : (length(y) == 2 ? sum(connects)^2 : Ah))
    end
end

function timing(operation::String, settings::Dict{String,Any})
    timing = string(operation," lasted: ")
    open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"a") do sim
        println(sim,timing)
    end
    #print to terminal
    if settings["timemsg"]
        println(timing)
    end
end

"""
simulate!()
"""
function simulate()
    # INITIALIZATION
    # Read in command line arguments
    settings = parse_commandline()

    srand(settings["rseed"])

    # unity test
    println(keys(settings))

    # Store landscape configuration
    landpars = read_landpars(settings)
    connects = readdlm(settings["connects"])
    # unity test
    println("Land init stored in object of type $(typeof(landpars))")

    # Store species information
    orgsref = read_spinput(settings)
    # unity test
    println("Sp info stored in object of type $(typeof(orgsref))")

    # Set insects implicit simulation
    interaction, scen, td, regime = implicit_insect(settings)
    
    # Create landscape
    mylandscape, landavail = landscape_init(landpars)
    # unity test
    println("Landscape initialized: type $(typeof(mylandscape))")

    # Initialize individual tagger: It tags all individuals ever created. It does not change if individuals die, so there is no risk of re-use of a tag.
    id_counter = 0
    
    # Create initial individuals
    orgs, id_counter = initorgs(landavail, orgsref, id_counter)
    
    println("Plants initialized: type $(typeof(orgs))")

    println("Starting simulation")

    try
        mkpath("$(settings["outputat"])/$(settings["simID"])")
        println("Output will be written to '$(settings["outputat"])'")
    catch
        println("Overwriting results to existing '$(settings["outputat"])/$(settings["simID"])' folder")
    end

    cd(pwd())

    # OUTPUT SIMULATION SETTINGS
    open(abspath(joinpath(settings["outputat"],settings["simID"],"simID")),"w") do ID
        println(ID, settings)
    end
    # START ID SIMULATION LOG FILE
    open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"w") do sim
        println(sim,string("Simulation: ",settings["simID"],now()))
    end
    # START SEED PRODUCTION FILE
    open(abspath(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv")),"w") do seedfile
        writedlm(seedfile, hcat(["week" "sp" "stage" "mode"], "abundance"))
    end
    
    # MODEL RUN
    for t in 1:settings["timesteps"]

        println("running week $t")
        
        # track simulation: header
        open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"a") do sim
            println(sim, "WEEK $t")
        end

        # UPDATE LANDSCAPE: weekly temperature and precipitation
        T = updateenv!(mylandscape, t, landpars)
        
        # DISTURBANCE
        if settings["disturb"] != "none"  
            tdist = disturb!(mylandscape,landavail,orgs,t,settings)
            # Initialize or update landscape properties that are relevant for life cycle processes
            if t == 1 || (settings["disturb"] in ["loss" "frag"] && t in [(tdist-1) tdist (tdist+1)])
                losschange(landavail, settings, t, tdist)
                fragchange(landavail, settings, t, connects, tdist)
            end
        end
        # OUTPUT: First thing, to see how community is initialized
        tic()
        orgstable(orgsref,landpars,orgs,mylandscape,t,settings)
        timing("WRITING ORGSOUTPUT", settings)
        toc()

        # LIFE CYCLE
        tic()
        if t == 1
            nogrowth = Int64[]
        end

        survive!(orgs,t,cK,K,settings,orgsref,landavail,T,nogrowth)

        global nogrowth = allocate!(mylandscape,orgs,t,aE,Boltz,settings,orgsref,T)

        develop!(orgs,orgsref)

        mate!(orgs,t,settings,scen,td,regime, 1)

        id_counter = mkoffspring!(orgs,t,settings,orgsref,id_counter)

        seedsi = release!(orgs,t,settings,orgsref)

        disperse!(landavail,seedsi,orgs,t,settings,orgsref,connects,AT,Ah)

        establish!(mylandscape,orgs,t,settings,orgsref,T)
        
        shedd!(orgs,orgsref,t)
        
        timing("LIFE CYCLE", settings)
        toc()
    end
end
simulate()
