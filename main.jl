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
using Organisms
using Setworld
using RCall

#const Boltz = 1.38064852e-23 # Alternatively:
const Boltz = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
#const aE = 1e-19
const aE = 0.63
global K = 0.0
global cK = 0.0
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

        "--initialland"
        help = "Name of file with landscape size values: areas of fragments, mean (and s.d.) temperature."
        arg_type = String
        default = abspath(pwd(),"inputs/initialland.jl") #"inputs/initialland.shp"

        "--disturbtype"
        help = "Type of environmental disturbance to be implemented: habitat area loss \"loss\", habitat fragmentation \"frag\" or temperature change \"temp\""
        arg_type = String
        default = "none"
        
        "--landmode"
        help = "Choose between using shape files to simulate the landscape and its change (\"real\") or providing the dimensions of the landscape to be simualted (total area, habitat area, number of habitat patches and distances between patches - \"artificial\")."
        arg_type = String
        default = "artif"
        
        "--landbuffer"
        help = "Buffer shape file or file containing its area"
        arg_type = String
        default = abspath(pwd(),"inputs/landbuffer.jl") # "inputs/landbuffer.shp"
        
        "--disturbland"
        help = "Either a shape file (if \`landmode\`)"
        arg_type = Any
        default = nothing #abspath(pwd(),"inputs/disturbland.jl") #abspath(pwd(),"inputs/disturbland.shp")

        "--tdist"
        help = "Timestep(s) where habitat loss or fragmentation is implemented."
        arg_type = String # path file to a vector of disturbance times
        
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
        Two methods because "real" landscapes do not need 
        """

function read_landpars(settings::Dict{String,Any})
    # Read in temperature time series, which is used in both modes
    temp_tsinput = loadtable(settings["temp_ts"])
    
    # Landscape is built differently, depending on the mode of simulation (which requires different types of files)

    if settings["landmode"] == "real"
        
        # send file names to R
        initialland = settings["initialland"]
        @rput initialland
        disturbland = settings["disturbland"]
        @rput disturbland
        landbuffer = settings["landbuffer"]
        @rput landbuffer
        
        # get patches/fragments areas and distances from shape/raster files
        landconfig = rcopy(Array{Any}, R"landconfig.R")

        # store them in landpars
        # no need to check for disturbance type because in the absence of disturbance-related files, R returns Nullable values
        landpars = Setworld.LandPars(landconfig[[1]],
                                     Fileprep.areatocell(landconfig[[2]]),
                                     landconfig[[3]],
                                     landconfig[[4]],
                                     landconfig[[5]],
                                     Fileprep.areatocell(landconfig[[6]]),
                                     landconfig[[7]],
                                     landconfig[[8]],
                                     select(temp_tsinput,:meantemp))

    elseif settings["landmode"] == "artif"

# deal with simulations where disturbance raster is needed
if settings["disturbtype"] in ["frag" "loss"]

        # send file names to R
        initialland = settings["initialland"]
        @rput initialland
        disturbland = settings["disturbland"]
        @rput disturbland
        
        # get patches/fragments areas and distances from raster files
        R"source(\"landnlmconfig.R\")"
        @rget initialmatrix
	@rget disturbmatrix

        landpars = Setworld.NeutralLandPars(Int.(initialmatrix),
        ifelse(disturbmatrix == "none", nothing, Int.(disturbmatrix)),
	select(temp_tsinput,:meantemp))
else
# send file names to R
        initialland = settings["initialland"]
        @rput initialland
        disturbland = settings["disturbland"]
        @rput disturbland
        
        # get patches/fragments areas and distances from raster files
        landconfig = rcopy(Array{Any}, R"source(\"landnlmconfig.R\")")
        @rget initialmatrix

        landpars = Setworld.NeutralLandPars(Int.(initialmatrix),
        nothing,
	select(temp_tsinput,:meantemp))
end
    else
        error("Please choose a mode of landscape simulation.")
    end
    
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
    td = select(insectsinput, :td) # time of perturbation
    regime = select(insectsinput, :regime)[1] # regime of loss
    remaining = select(insectsinput, :remaining) # proportion of pollination service loss
    
    return interaction, scen, td, regime, remaining
end

"""
        outputorgs(orgs,t,settings)
        Saves a long format table with the organisms field informations.
        """
function orgstable(orgs::Array{Organisms.Organism,1}, t::Int64, settings::Dict{String,Any})

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
        
        # TODO: output landscape 
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
        
        #TODO output landscape
    end
end

"""
        loaddisturbance()
        Store parameters necessary to implement disturbance.
        """
function loaddisturbance(settings)

tdist = 0
    # read in the disturbance file, if not done so yet

    # select file according to keyword: loss, frag, temp
    if settings["disturbtype"] != "none"
        if settings["disturbtype"] in ["loss" "frag"]
            
            if settings["landmode"] == "real"
                tdist = select(loadtable(settings["tdist"]), :tdist)
                #disturblandspecs = nothing
            elseif settings["landmode"] == "artif"
                # this mode also requires a file specifing the paths to the connectivity matrices
                tdist = select(loadtable(settings["tdist"]), :tdist)
                #include(settings["disturbland"]) # unnecessary for continuous landscapes because it will be a file path to the raster file, which is read in read_landpars()
                #disturblandspecs = [fragsid, fareas, disturbconnect]
            end


elseif settings["disturbtype"] == "poll"
	    tdist = select(loadtable(settings["insect"]), :td)
            println("Pollination loss is simulated according to parameters in the \'insects\' file.")
            
        elseif settings["disturbtype"] == "temp"
            println("Temperature change is simulated with the temperature file provided.")
        else
            error("Please specify one of the disturbance scenarios with `--disturb`:
                          \n\'none\' if no disturbance should be simulated,
                          \n\'loss\' for habitat area loss,
                          \n\'frag\' for habitat fragmentation,
                          \n\'temp\' for temperature change,
                          \n\'poll\' for pollination loss.")
        end
    end

return tdist
end

"""
        disturb!()

        """
function disturb!(landscape::Array{Dict{String,Float64},N} where N, landavail::BitArray{N} where N, orgs::Array{Organisms.Organism,1}, t::Int64, tdist::Array{Int64,1}, settings::Dict{String,Any}, landpars::NeutralLandPars)
    
        if settings["disturbtype"] == "loss"
            loss = prod(size(landpars.disturbland))/prod(landpars.initialarea)
            landavail = Setworld.destroyarea!(landpars, landavail, settings)
            Organisms.destroyorgs!(orgs, landavail, settings)
        elseif settings["disturbtype"] == "frag"
            landscape, landavail = Setworld.disturbland!(landscape, landavail, landpars)
            Organisms.destroyorgs!(orgs, landavail, settings)
        end                

return landscape,landavail

end
	
"""
        updateK!()
        Updates the carrying capacity of the landscape (`K`) and of each gridcell (`cK`). 
        """

function updateK!(landavail::BitArray{2}, settings::Dict{String,Any}, t::Int64, tdist::Any)

    # check on which timesteps to write land dims (ifelse() does not work)
    if tdist == nothing
       criticalts = [1]
    else
	criticalts = [1 (tdist-1) tdist (tdist+1)]
    end
    
    if t in criticalts
        # output message
        if t == 1
            open(abspath(joinpath(settings["outputat"],settings["simID"],"landlog.txt")),"w") do sim
                println(sim, "week\tK\tcK\thabitatarea\ttotalarea")
            end          
        end

        # habitat area
        habitatarea = length(find(x -> x == true, landavail))*625 # cells area 25x25cm
	totalarea = prod(size(landavail))
	
        # K and grid-cell K
        global K = (1.5/100)*habitatarea # xtons/ha = x.10-²g/1m²
        global cK = K/length(find(x -> x == true, landavail))

        open(abspath(joinpath(settings["outputat"],settings["simID"],"landlog.txt")),"a") do sim
            writedlm(sim,[t K cK habitatarea totalarea])
        end
    end
    
end

function timing(operation::String, settings::Dict{String,Any})
    timing_stamp = string(operation," lasted: ")
    open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"a") do sim
        println(sim,timing_stamp)
    end
    #print to terminal
    if settings["timemsg"]
        println(timing_stamp)
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

    # Load disturbance parameters, if necessary
    if settings["disturbtype"] != "none"
        tdist = loaddisturbance(settings)
    else
        tdist = nothing 
    end            
    
    # Store landscape configuration
    landpars = read_landpars(settings)
    
    # unity test
    println("Land init stored in object of type $(typeof(landpars))")

    # Store species information
    orgsref = read_spinput(settings)
    # unity test
    println("Sp info stored in object of type $(typeof(orgsref))")

    # Set insects implicit simulation
    interaction, scen, td, regime, remaining = implicit_insect(settings)
    
    # Create landscape
    mylandscape, landavail = Setworld.landscape_init(landpars)
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

    ##############################
    # OUTPUT SIMULATION SETTINGS #
    ##############################
    open(abspath(joinpath(settings["outputat"],settings["simID"],"simsettings.jl")),"w") do ID
        println(ID, "tdist = ", tdist)
	println(ID, "landpars = ", typeof(landpars), "\ninitial = ", typeof(landpars.initialland), "\ndisturb = ", typeof(landpars.disturbland))  
	println(ID, "interaction = ", interaction)
	println(ID, "scen = ", scen)
	println(ID, "td = ", td)
	println(ID, "regime = ", regime)
	println(ID, "remaining = ", remaining)
	println(ID, "commandsettings = ", settings)
    end
    
    # START ID SIMULATION LOG FILE
    open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"w") do sim
        println(sim,string("Simulation: ",settings["simID"],now()))
    end

    # START SEED PRODUCTION FILE
    open(abspath(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv")),"w") do seedfile
        writedlm(seedfile, hcat(["week" "sp" "stage" "mode"], "abundance"))
    end

    #############
    # MODEL RUN #
    #############
    for t in 1:settings["timesteps"]

        println("running week $t")
        
        # track simulation: header
        open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"a") do sim
            println(sim, "WEEK $t")
        end

        # UPDATE LANDSCAPE & CARRYING CAPACITY: weekly temperature and precipitation
        T = updateenv!(t, landpars)
        updateK!(landavail, settings, t, tdist)
	
        # IMPLEMENT LANDSCAPE DISTURBANCE
        if settings["disturbtype"] in ["frag" "loss"] && t in tdist  
            landscape, landavail = disturb!(mylandscape,landavail,orgs,t,tdist,settings,landpars)
        end
        
        # OUTPUT: First thing, to see how community is initialized
        tic()
        orgstable(orgs,t,settings)
        timing("WRITING ORGSOUTPUT", settings)
        toc()

        # LIFE CYCLE
        tic()

        survive!(orgs, t, cK, K, settings, orgsref, landavail, T, nogrowth)

        global nogrowth = allocate!(orgs, t, aE, Boltz, settings, orgsref, T)

        develop!(orgs, orgsref)

        mate!(orgs, t, settings, scen, regime, td, remaining)

        id_counter = mkoffspring!(orgs, t, settings, orgsref, id_counter, landavail)

        seedsi = release!(orgs, t, settings, orgsref) # only recently released seeds need to disperse. The others only need to survive

        disperse!(landavail, seedsi, orgs, t, settings, orgsref, landpars, tdist)

        establish!(orgs, t, settings, orgsref, T)
        
        shedd!(orgs, orgsref, t)
        
        timing("LIFE CYCLE", settings)
        toc()
    end
end
simulate()
