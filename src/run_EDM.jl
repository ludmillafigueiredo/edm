#!/usr/bin/env julia

# Get model directory and include it in Julia's loading path
#
#cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/") # Only in ATOM
EDDir = abspath(joinpath(pwd(), "src"))
push!(LOAD_PATH,EDDir)

# Load Julia & model packages
using ArgParse
using Distributions
using JuliaDB #for in/outputs
using CSV #TODO replace JuliaDB with it?
using DataFrames
using DataValues
using RCall

# Load EDM modules
using auxfunctions
using submodels
using outputs

#const Boltz = 1.38064852e-23 # Alternatively:
const Boltz = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
#const aE = 1e-19
const aE = 0.63
global K = 0.0
global cK = 0.0

function parse_commandline()
    sets = ArgParseSettings() #object that will be populated with the arguments by the macro

    @add_arg_table sets begin
        "--simID"
        help = "Name of the folder where outputs will be stored."
        arg_type = String
        required = true

        "--nreps"
        help = "Name of the folder where outputs will be stored."
        arg_type = Int64
        default = 1

        "--rseed"
        help = "Seed for RNG"
        arg_type = Int
        required = true

        "--outputat"
        help = "Name of directory where output should be written ."
        arg_type = String
        default = abspath(pwd(),"outputs")

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
        help = "Output timing to terminal, as well as checkpoint"
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
        disturbland = CSV.read(settings["disturbland"], header = true, types = Dict("td" => Int64, "disturbland" => Any))[:disturbland]
        @rput disturbland
        landbuffer = settings["landbuffer"]
        @rput landbuffer

        # get patches/fragments areas and distances from shape/raster files
        landconfig = rcopy(Array{Any}, reval(abspath(joinpath(EDDir,"landconfig.R"))))

        # store them in landpars
        # no need to check for disturbance type because in the absence of disturbance-related files, R returns Nullable values
        landpars = submodels.LandPars(landconfig[[1]],
                                     auxfunctions.areatocell(landconfig[[2]]),
                                     landconfig[[3]],
                                     landconfig[[4]],
                                     landconfig[[5]],
                                     auxfunctions.areatocell(landconfig[[6]]),
                                     landconfig[[7]],
                                     landconfig[[8]],
                                     select(temp_tsinput,:meantemp))

    elseif settings["landmode"] == "artif"

        # deal with simulations where disturbance raster is needed
        if settings["disturbtype"] in ["frag" "loss"]

            # send file names to R
            initialland = settings["initialland"]
            @rput initialland
            # assigning type Any to disturbland is not possible
            disturbland = settings["disturbland"] #object has to be initialized
            if settings["disturbtype"] == "loss"
                disturbland = CSV.read(settings["disturbland"], header = true, types = Dict("td" => Int64, "proportion" => Float64))
            elseif settings["disturbtype"] == "frag"
                disturbland = CSV.read(settings["disturbland"], header = true, types = Dict("disturbland" => String))
            end
            @rput disturbland

            # get 0-1 habitat suitability matrices
            landconfig = rcopy(Array{Any}, reval(string("source(\"",abspath(joinpath(EDDir,"landconfig.R")),"\")")))
            @rget initialmatrix
            @rget disturbmatrix

            # if fragmentation is simulated, the file names are assgined to dist field, if area loss, the proportion holding columns, if none, default value nothing
            if disturbmatrix == "notfrag"
                landpars = submodels.NeutralLandPars(Int.(initialmatrix),
                                                    disturbland,
                                                    select(temp_tsinput,:meantemp))
            else
                landpars = submodels.NeutralLandPars(Int.(initialmatrix),
                                                    Int.(disturbmatrix),
                                                    select(temp_tsinput,:meantemp))
            end
            #TODO block above can handle pollination
        else
            # send file names to R
            initialland = settings["initialland"]
            @rput initialland
            disturbland = settings["disturbland"]
            @rput disturbland

            # get suitability matrices from R
            landconfig = rcopy(Array{Any}, reval(string("source(\"",abspath(joinpath(EDDir,"landnlmconfig.R")),"\")")))
            @rget initialmatrix
            @rget disturbmatrix

            landpars = submodels.NeutralLandPars(Int.(initialmatrix),
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
                Reads in species initial conditions and parameters. Stores tehm in `sppref`, a structure with parameters names as Dictionnary fields, where species names are the keys to the parameter values.
                """

function read_spinput(settings::Dict{String,Any})

    #spinputtbl = loadtable(abspath(pwd(),"inputs/species.csv"))
    spinputtbl = loadtable(settings["spinput"])

    sppref = SppRef(Array(rows(spinputtbl,:sp_id)),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:abund)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:kernel)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:clonality)[i] == "true" #translates R's true into Julia's TRUE
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:seedmass)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:maxmass)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:span_min)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:span_max)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:firstflower_min)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:firstflower_max)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:floron)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:floroff)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:seednumber_min)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:seednumber_max)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:seedon)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:seedoff)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:bankduration_min)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:bankduration_max)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:b0grow)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:b0germ)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:b0mort)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:temp_opt)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict(rows(spinputtbl,:sp_id)[i] =>
                                    rows(spinputtbl,:temp_tol)[i]
                                for i in 1:length(rows(spinputtbl,:sp_id))),
                           Dict())

    return sppref
end

function define_traitranges(settings::Dict{String,Any})

    #spinputtbl = loadtable(abspath(pwd(),"inputs/species.csv"))
    spinputtbl = loadtable(settings["spinput"])

    traitranges = TraitRanges(
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [0.95*rows(spinputtbl,:seedmass)[i], 1.05*rows(spinputtbl,:seedmass)[i]]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [0.95*rows(spinputtbl,:maxmass)[i], 1.05*rows(spinputtbl,:maxmass)[i]]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(0.95*rows(spinputtbl,:span_min)[i], RoundDown)), Int(round(1.05*rows(spinputtbl,:span_max)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(0.95*rows(spinputtbl,:firstflower_min)[i], RoundDown)), Int(round(1.05*rows(spinputtbl,:firstflower_max)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(0.95*rows(spinputtbl,:floron)[i], RoundDown)), Int(round(1.05*rows(spinputtbl,:floron)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(0.95*rows(spinputtbl,:floroff)[i], RoundDown)), Int(round(1.05*rows(spinputtbl,:floroff)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(0.95*rows(spinputtbl,:seednumber_min)[i], RoundDown)), Int(round(1.05*rows(spinputtbl,:seednumber_max)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(0.95*rows(spinputtbl,:seedon)[i], RoundDown)), Int(round(1.05*rows(spinputtbl,:seedon)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(0.95*rows(spinputtbl,:seedoff)[i], RoundDown)), Int(round(1.05*rows(spinputtbl,:seedoff)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(0.95*rows(spinputtbl,:bankduration_min)[i], RoundDown)), Int(round(1.05*rows(spinputtbl,:bankduration_max)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id)))
    )

    return traitranges

end

"""
                implicit_insect(settings)
                Reads how insects are going to be implicitly simulated.
                """
function implicit_insect(settings::Dict{String,Any})

    insectsinput = loadtable(settings["insect"])

    interaction = select(insectsinput, :interaction)[1]
    scen = select(insectsinput, :scen)[1] # pollination scenario
    remaining = select(insectsinput, :remaining) # proportion of pollination service loss

    return interaction, scen, remaining
end

"""
                outputorgs(plants,t,settings)
                Saves a long format table with the organisms field informations.
                """
function orgstable(plants::Array{submodels.Plant,1}, t::Int64, settings::Dict{String,Any})

    outplants = find(x -> (x.stage in ["j" "a"] || (x.stage == "s" && x.age > 1)), plants)

    # output header
    if t == 1
        header = hcat(["week"],
                      reshape(string.(fieldnames(Plant)[1:20]),1,:),
                      ["leaves" "stem" "root" "repr"],
                      reshape(string.(fieldnames(Plant)[22:end]),1,:))
        open(abspath(joinpath(settings["outputat"],settings["simID"],"orgsweekly.txt")), "w") do output
            writedlm(output, header) #reshape(header, 1, length(header)))
        end
    end

    # output plants info
    if t == 1 || rem(t,settings["tout"]) == 0

        for o in outplants
            open(abspath(joinpath(settings["outputat"],settings["simID"],"orgsweekly.txt")), "a") do output
                writedlm(output, hcat(t,
                                      plants[o].id,
                                      plants[o].stage,
                                      plants[o].location,
                                      plants[o].sp,
                                      plants[o].kernel,
                                      plants[o].clonality,
                                      plants[o].seedmass,
                                      plants[o].maxmass,
                                      plants[o].span,
                                      plants[o].firstflower,
                                      plants[o].floron,
                                      plants[o].floroff,
                                      plants[o].seednumber,
                                      plants[o].seedon,
                                      plants[o].seedoff,
                                      plants[o].bankduration,
                                      plants[o].b0grow,
                                      plants[o].b0germ,
                                      plants[o].b0mort,
                                      plants[o].age,
                                      plants[o].mass["leaves"],
                                      plants[o].mass["stem"],
                                      plants[o].mass["root"],
                                      plants[o].mass["repr"],
                                      plants[o].mated))
            end
        end
    end
end

"""
                loaddisturbance()
                Store parameters necessary to implement disturbance.
                """
function loaddisturbance(settings)

    tdist = nothing
    # read in the disturbance file, if not done so yet

    # select file according to keyword: loss, frag, temp
    if settings["disturbtype"] != "none"
        if settings["disturbtype"] in ["loss" "frag"]

            tdist = CSV.read(settings["disturbland"], header = true, types = Dict("td" => Int64))[:td]

            return tdist

        elseif settings["disturbtype"] == "poll"
            tdist = select(loadtable(settings["insect"]), :td)
            println("Pollination loss is simulated according to parameters in the \'insects\' file.")
            return tdist

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
end

"""
                disturb!()

                """
function disturb!(landscape::Array{Dict{String,Float64},N} where N, landavail::BitArray{2}, plants::Array{submodels.Plant,1}, t::Int64, settings::Dict{String,Any}, landpars::NeutralLandPars, tdist::Any)

    if settings["disturbtype"] == "loss"
        landscape, landavail = submodels.destroyarea!(landpars, landavail, settings, t)
    elseif settings["disturbtype"] == "frag"
        landscape, landavail = submodels.fragment!(landscape, landavail, landpars, t, tdist)
    end

    submodels.destroyorgs!(plants, landavail, settings)

    return landscape,landavail

end

"""
                updateK!()
                Updates the carrying capacity of the landscape (`K`) and of each gridcell (`cK`).
                """

function updateK!(landavail::BitArray{2}, settings::Dict{String,Any}, t::Int64, tdist::Any)

    criticalts = Array{Int64,N} where N
    # check on which timesteps to write land dims (ifelse() does not work)
    if tdist == nothing
        criticalts = [1]
    else
        criticalts = vcat(1, (tdist-1), tdist, (tdist+1))
    end

    if t in criticalts
        # output message
        if t == 1
            open(abspath(joinpath(settings["outputat"],settings["simID"],"landlog.txt")),"w") do sim
                println(sim, "week\tK\tcK\thabitatarea\ttotalarea")
            end
        end

        # habitat area
        habitatarea = length(find(x -> x == true, landavail))*10000 # cells area 100x100cm; habitatarea in cm²
        totalarea = prod(size(landavail))

        # K and grid-cell K
        global K = (3.5/100)*habitatarea # xtons/ha = x.10-2g/1cm²
        global cK = (3.5/100)*10000	 

        open(abspath(joinpath(settings["outputat"],settings["simID"],"landlog.txt")),"a") do sim
            writedlm(sim,[t K cK habitatarea totalarea])
        end
    end

end

function updateK!(landavail::BitArray{2}, settings::Dict{String,Any}, t::Int64)

    # habitat area

    habitatarea = length(find(x -> x == true, landavail))*10000 # cells area 25x25cm; habitatarea in cm²
    totalarea = prod(size(landavail))

    # K and grid-cell K
    global K = (3.5/100)*habitatarea # xtons/ha = x.10-2g/1cm²
    globalcK = K/length(find(x -> x == true, landavail))
end


function timing(operation::String, settings::Dict{String,Any})
    timing_stamp = string(operation," lasted: ")
    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        println(sim,timing_stamp)
    end
    #print to terminal
    if settings["timemsg"]
        println(timing_stamp)
    end
end

"""
            updatefitness!(sppref, mean_opt::Float64, std_tol::Float64, mean_annual::Float64, max_fitness::Float64)

            Calculate the fitness value according to a Gauss function:

                f(x) = a*exp(-((x-b)²)/(2*c²)),

            # Arguments
            - `fitnessdict::Dict{String, Array{Float64, 1}}`: dictionnary holding the species performance for each timestep
            - `maximal fitness::Float64=1.0`: parameter `a`, the height of the curve's peak.
            - `mean_opt::Float64`: parameter `b` is the position of the center of the peak.
            - `std_tol::Float64`: parameter `c` is the standard deviation, 
             """

function updatefitness!(sppref::SppRef, mean_annual::Float64, max_fitness::Float64)

    for sp in sppref.sp_id
        
        mean_opt = sppref.temp_opt[sp]
        std_tol = sppref.temp_tol[sp]

        absolute_fitness = max_fitness*exp(-((mean_annual-mean_opt)^2)/(2*(std_tol^2)))
        sppref.fitness[sp] = absolute_fitness
    end

end


function updatefitness!(sppref::SppRef, mean_annual::Float64, max_fitness::Float64, t::Int64, settings::Dict{String, Any})

    for sp in sppref.sp_id
        
        mean_opt = sppref.temp_opt[sp]
        std_tol = sppref.temp_tol[sp]

        absolute_fitness = max_fitness*exp(-((mean_annual-mean_opt)^2)/(2*(std_tol^2)))
        sppref.fitness[sp] = absolute_fitness
    end

    open(abspath(joinpath(settings["outputat"], settings["simID"], "spp_fitness.csv")),"a") do fitnessfile
    	writedlm(fitnessfile, hcat(t, sp, get(sppref.fitness, sp, "NA")) for sp in collect(keys(sppref.fitness)))
    end
end


"""
                simulate!()
                """
function simulate()

    # INITIALIZATION 
    #################

    # Read in command line arguments
    settings = parse_commandline()

    # unity test
    println(keys(settings))

    # Load disturbance parameters, if necessary
    tdist = nothing
    if settings["disturbtype"] != "none"
        tdist = loaddisturbance(settings)
    end
    println(tdist)

    # Store landscape configuration
    landpars = read_landpars(settings)

    # unity test
    println("Land init stored in object of type $(typeof(landpars))")

    # Store species information
    sppref = read_spinput(settings)
    ## unity test
    println("Sp info stored in object of type $(typeof(sppref))")
    traitranges = define_traitranges(settings)
    ## unity test
    println("Sp ranges stored in object of type $(typeof(traitranges))")

    # Set insects implicit simulation
    interaction, scen, remaining = implicit_insect(settings)

    # Create landscape
    mylandscape, landavail = submodels.landscape_init(landpars)
    # unity test
    println("Landscape initialized: type $(typeof(mylandscape))")

    # Initialize individual tagger: It tags all individuals ever created. It does not change if individuals die, so there is no risk of re-use of a tag.
    id_counter = 0

    # Initialize management counter
    management_counter = 0

    # Initialize abundances according to species fitness
    updateK!(landavail, settings, 1)
    T, mean_annual = updateenv!(1, landpars)
    updatefitness!(sppref, mean_annual, 1.0)

    # Create initial individuals
    plants, id_counter = initorgs(landavail, sppref, id_counter, settings, K)

    println("Plants initialized: type $(typeof(plants))")

    cd(pwd())

    println("Starting simulation")

    # ORGANIZE OUTPUT FOLDERS
    #########################
    
    results_folder = string()
    for rep in 1:settings["nreps"]

        srand(settings["rseed"])

        if settings["nreps"] == 1
            simresults_folder = joinpath(settings["outputat"], settings["simID"])
            results_folder = simresults_folder
        else
            simresults_folder = joinpath(settings["outputat"], string(settings["simID"], "_", rep))
            results_folder = joinpath(settings["outputat"], settings["simID"])
        end

        try
            mkpath("$(simresults_folder)")
            println("Output will be written to '$(simresults_folder)'")
        catch
            println("Overwriting results to existing '$(simresults_folder)' folder")
        end

        # OUTPUT SIMULATION SETTINGS 
        open(abspath(joinpath(simresults_folder, "simsettings.jl")),"w") do ID
            println(ID, "tdist = ", tdist)
            println(ID, "landpars = ", typeof(landpars), "\ninitial = ", typeof(landpars.initialland), "\ndisturb = ", typeof(landpars.disturbland))
            println(ID, "interaction = ", interaction)
            println(ID, "scen = ", scen)
            println(ID, "remaining = ", remaining)
            println(ID, "commandsettings = ", settings)
        end

        # INITIALIZE FILE TO LOG SIMULATION PROGRESS
        open(abspath(joinpath(simresults_folder, "checkpoint.txt")),"w") do sim
            println(sim,string("Simulation: ",settings["simID"],now()))
        end

        # INITIALIZE FILE TO LOG LIFE-HISTORY EVENTS
        open(abspath(joinpath(simresults_folder, "eventslog.txt")),"w") do sim
            writedlm(sim, hcat("week", "event", "stage", "age"))
        end

        # INITIALIZE FILE TO LOG METABOLIC RATES
        open(abspath(joinpath(settings["outputat"],settings["simID"],"metaboliclog.txt")),"a") do sim

            writedlm(sim, hcat("stage", "age", "rate", "probability", "event"))

        end

        # INITIALIZE FILE TO LOG SEED PRODUCTION
        open(abspath(simresults_folder, "offspringproduction.csv"),"w") do seedfile
            writedlm(seedfile, hcat(["week" "sp" "stage" "mode"], "abundance"))
        end

	# INITIALIZE FILE TO LOG SPECIES FITNESS
        open(abspath(simresults_folder, "spp_fitness.csv"),"w") do fitnessfile
            writedlm(fitnessfile, hcat("week", "sp", "fitness"))
        end

        # MODEL RUN 
        ############
        
        for t in 1:settings["timesteps"]

            # check-point
            open(abspath(joinpath(simresults_folder, "checkpoint.txt")),"a") do sim
                println(sim, "\nWEEK $t")
		println(sim, "Species richness: $(length(unique(map(x -> x.sp, plants))))")
            end
	    
	    println("\nWEEK $t")
	    println("Species richness: $(length(unique(map(x -> x.sp, plants))))")

            # IMPLEMENT LANDSCAPE DISTURBANCE
            if settings["disturbtype"] in ["frag" "loss"] && t in tdist
                landscape, landavail = disturb!(mylandscape,landavail,plants,t,settings,landpars,tdist)
            end
            updateK!(landavail, settings, t, tdist)

	    # UPDATE CURRENT TEMPERATURE
	    if rem(t, 52) == 1
                T, mean_annual = updateenv!(t, landpars)
	    else
	        T = updateenv!(t, landpars)
	    end
	    
            # UPDATE SPECIES FITNESS
            updatefitness!(sppref, mean_annual, 1.0, t, settings)
            
            # OUTPUT: First thing, to see how community is initialized
            tic()
            orgstable(plants,t,settings)
            timing("Time writing output", settings)
            toc()

            # MANAGEMENT: it happens at least once every year, annually
            if ((31 < rem(t,52) < 39) && management_counter < 1) #mowing cannot happen before the 1st of July
                management_counter = manage!(plants, t, management_counter, settings)
            end

            # UPDATE CURRENT BIOMASS PRODUCTION
            biomass_production = sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants), 0.00001))
            # check-point
            open(abspath(joinpath(simresults_folder, "checkpoint.txt")),"a") do sim
                writedlm(sim, hcat("Biomass production:", biomass_production))
            end

            # LIFE CYCLE
            tic()

            allocate!(plants, t, aE, Boltz, settings, sppref, T, biomass_production, K, "a")
            
            survive!(plants, t, cK, K, settings, sppref, landavail, T, biomass_production, "a")

	    allocate!(plants, t, aE, Boltz, settings, sppref, T, biomass_production, K, "j")
            
            survive!(plants, t, cK, K, settings, sppref, landavail, T, biomass_production, "j")

            develop!(plants, sppref, settings, t)

            mate!(plants, t, settings, scen, tdist, remaining)

            id_counter = mkoffspring!(plants, t, settings, sppref, id_counter, landavail, T, traitranges)

            seedsi = release!(plants, t, settings, sppref) # only recently released seeds need to disperse. The others only need to survive

            disperse!(landavail, seedsi, plants, t, settings, sppref, landpars, tdist)

            establish!(plants, t, settings, sppref, T, biomass_production, K)
            
            shedd!(plants, sppref, t, settings)

            timing("Time running life cycle:", settings)
            toc()
        end
end

return settings, results_folder

end

# run simulation
settings, results_folder = simulate()

# analyse results
analysED(settings, results_folder)
