#!/usr/bin/env julia

EDDir = joinpath(pwd(), "src") # TODO test without pwd()
push!(LOAD_PATH,EDDir)

using ArgParse
using Distributions
using JuliaDB
using CSV #TODO replace JuliaDB with it?
using DataFrames
using DataValues
using RCall
using auxfunctions
using entities
using submodels
using outputs

include("constants_globalpars.jl")

function parse_commandline()
    sets = ArgParseSettings() #object that will be populated with the arguments by the macro
    @add_arg_table sets begin
        "--simID"
        help = "Name of the folder where outputs will be stored."
        arg_type = String
        required = true

        "--nreps"
        help = "Number of replicates."
        arg_type = Int64
        default = 1

        "--rseed"
        help = "Seed for RNG"
        arg_type = Int
        required = true

        "--outputat"
        help = "Name of directory where output should be written ."
        arg_type = String
        default = "outputs"

        "--spinput"
        help = "Name of file with species list."
        arg_type = String
        default = joinpath("inputs", "species.csv")

        "--insect"
        help = "How to explicitly model insects:
                pollination-independent reproduction \"indep\";
                equal pollination loss for all species \"equal\"."
        arg_type = String
        default = joinpath("inputs", "insects.csv")

        "--landbuffer"
        help = "Buffer shape file or file containing its area"
        arg_type = String
        default = joinpath("inputs", "landbuffer.jl")

        "--initialland"
        help = "Name of file with landscape size values: areas of fragments, mean (and s.d.) temperature."
        arg_type = String
        default = joinpath("inputs","initialland.jl")

        "--disturbtype"
        help = "Type of environmental disturbance to be implemented: habitat area loss \"loss\", habitat fragmentation \"frag\" or temperature change \"temp\""
        arg_type = String
        default = "none"

        "--landmode"
        help = "Choose between using shape files to simulate the landscape and its change (\"real\") or providing the dimensions of the landscape to be simualted (total area, habitat area, number of habitat patches and distances between patches - \"artificial\")."
        arg_type = String
        default = "artif"

        "--disturbland"
        help = "Either a shape file (if \`landmode\`)"
        arg_type = Any
        default = nothing

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
        default = joinpath("inputs", "temp1917_2017.csv")

        "--timemsg"
        help = "Output timing to terminal, as well as checkpoint"
        arg_type = Bool
        default = false
    end
    return parse_args(sets)
end


"""
    read_landin(settings)
Read, derive and store values related to the environmental conditions: temperature, landscape size, availability, and changes thereof.
"""
function read_landpars(settings::Dict{String,Any})
    # Temperature time-series
    temp_tsinput = loadtable(settings["temp_ts"])

    # Landscape is built differently, depending on the mode of simulation
    if settings["landmode"] == "real" # simulation arena is built based on raster files, which are analysed in R

        # send file paths to R
        initialland = settings["initialland"]
        @rput initialland
        # disturbland file contains time of habitat destruction and the raster file representing it afterwards
        disturbland = CSV.read(settings["disturbland"], header = true, types = Dict("td" => Int64, "disturbland" => Any))[:disturbland]
        @rput disturbland
        # buffer area
        landbuffer = settings["landbuffer"]
        @rput landbuffer

        # get patches/fragments areas and distances from shape/raster files
        landconfig = rcopy(Array{Any}, reval(joinpath(EDDir,"landconfig.R")))
        # store them in landpars
        # in the absence of disturbance-related files, R returns Nullable values
        landpars = submodels.LandPars(landconfig[[1]],
                                     auxfunctions.areatocell(landconfig[[2]]),
                                     landconfig[[3]],
                                     landconfig[[4]],
                                     landconfig[[5]],
                                     auxfunctions.areatocell(landconfig[[6]]),
                                     landconfig[[7]],
                                     landconfig[[8]],
                                     select(temp_tsinput,:meantemp))

    elseif settings["landmode"] == "artif" # simulation arena is built from a neutral landscape

        if settings["disturbtype"] in ["frag" "loss"]

            # send file names to R
            initialland = settings["initialland"]
            @rput initialland

            disturbland = settings["disturbland"] #object has to be initialized
            if settings["disturbtype"] == "loss"
                disturbland = CSV.read(settings["disturbland"], header = true, types = Dict("td" => Int64, "proportion" => Float64))
            elseif settings["disturbtype"] == "frag"
                disturbland = CSV.read(settings["disturbland"], header = true, types = Dict("disturbland" => String))
            end
            @rput disturbland

            # get 0-1 habitat suitability matrices
            landconfig = rcopy(Array{Any}, reval(string("source(\"",joinpath(EDDir,"landconfig.R"),"\")")))
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
            landconfig = rcopy(Array{Any}, reval(string("source(\"",joinpath(EDDir,"landnlmconfig.R"),"\")")))
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
Reads in species trait values (min. and max. values for the evolvable ones). Stores them in `sppref`, a structure with parameters as Dictionnary fields, where species names are the keys to the parameter values.
"""
function read_spinput(settings::Dict{String,Any})

    spinputtbl = loadtable(settings["spinput"])

    sppref = submodels.SppRef(Array(rows(spinputtbl,:sp_id)),
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
                                   rows(spinputtbl,:compartsize)[i]
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

"""
    define_traitranges(settings)
Based on input species traits values, set the limits upon which they can vary during sexual reproduction: [min., max.].
Values are store in an object of type TraitRanges, where each field is a Dictionary for a trait, where the keys are the species, and the values are Arrays holding min. and max. values the trait can take.
"""
function define_traitranges(settings::Dict{String,Any})

    spinputtbl = loadtable(settings["spinput"])

    traitranges = TraitRanges(
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [rows(spinputtbl,:seedmass)[i], rows(spinputtbl,:seedmass)[i]]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [rows(spinputtbl,:compartsize)[i], rows(spinputtbl,:compartsize)[i]]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(rows(spinputtbl,:span_min)[i], RoundDown)), Int(round(rows(spinputtbl,:span_max)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(rows(spinputtbl,:firstflower_min)[i], RoundDown)), Int(round(rows(spinputtbl,:firstflower_max)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(rows(spinputtbl,:floron)[i], RoundDown)), Int(round(rows(spinputtbl,:floron)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(rows(spinputtbl,:floroff)[i], RoundDown)), Int(round(rows(spinputtbl,:floroff)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(rows(spinputtbl,:seednumber_min)[i], RoundDown)), Int(round(rows(spinputtbl,:seednumber_max)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(rows(spinputtbl,:seedon)[i], RoundDown)), Int(round(rows(spinputtbl,:seedon)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(rows(spinputtbl,:seedoff)[i], RoundDown)), Int(round(rows(spinputtbl,:seedoff)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))),
        Dict(rows(spinputtbl,:sp_id)[i] =>
             [Int(round(rows(spinputtbl,:bankduration_min)[i], RoundDown)), Int(round(rows(spinputtbl,:bankduration_max)[i], RoundUp))]
             for i in 1:length(rows(spinputtbl,:sp_id))))
    return traitranges
end

"""
    implicit_insect(settings)
Reads how insects are going to be implicitly simulated.
"""
function implicit_insect(settings::Dict{String,Any})

    insectsinput = loadtable(settings["insect"])
    interaction = select(insectsinput, :interaction)[1] # type of interaction which is affected
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
        open(joinpath(settings["outputat"],settings["simID"],"orgsweekly.txt"), "w") do output
            writedlm(output, header) #reshape(header, 1, length(header)))
        end
    end

    # output plants info
    if t == 1 || rem(t,settings["tout"]) == 0

        for o in outplants
            open(joinpath(settings["outputat"],settings["simID"],"orgsweekly.txt"), "a") do output
                writedlm(output, hcat(t,
                                      plants[o].id,
                                      plants[o].stage,
                                      plants[o].location,
                                      plants[o].sp,
                                      plants[o].kernel,
                                      plants[o].clonality,
                                      plants[o].seedmass,
                                      plants[o].compartsize,
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
    disturb!(landscape, landavail, plants, t, settings, landpars, tdist)
Change landscape structure and metrics according to the simulation scenario (`loss`, habitat loss, or `frag`, fragmentation)
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
    initK(landavail, settings, t)
Calculate the carrying capacity of the landscape (`K`) and of each gridcell (`cK`) at initialization.
`K` is used to initialize the species with abundances corresponding to a niche partitioning model.
"""
function initK(landavail::BitArray{2}, settings::Dict{String,Any}, t::Int64)
    habitatarea = length(find(x -> x == true, landavail))
    totalarea = prod(size(landavail))
    K = cK*habitatarea
    return K
end

"""
    updateK!(landavail, settings, t, tdist)
Update the carrying capacity of the landscape (`K`) and of each gridcell (`cK`).
It is called during initialization (`t` = 1)  and is called again if landscape is disturbed (at `tdist`; `criticalts` keeps track of the updates and outputs the new values).
"""
function updateK!(K::Float64, landavail::BitArray{2}, settings::Dict{String,Any}, t::Int64, tdist::Any)

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
            open(joinpath(settings["outputat"],settings["simID"],"landlog.txt"),"w") do sim
                println(sim, "week\tK\tcK\thabitatarea\ttotalarea")
            end
        end

        # habitat area
        habitatarea = length(find(x -> x == true, landavail))  #number of habitat grid cels ggrid cells
        totalarea = prod(size(landavail)) # total number of grid cell. Habitat or not

	K = cK*habitatarea
        open(joinpath(settings["outputat"],settings["simID"],"landlog.txt"),"a") do sim
            writedlm(sim,[t K cK habitatarea totalarea])
        end
    end
    return K
end

"""
    updatefitness!(sppref, mean_opt::Float64, std_tol::Float64, mean_annual::Float64, max_fitness::Float64)

Calculate the fitness value according to a Gauss function:

         f(x) = a*exp(-((x-b)²)/(2*c²)),

and store the value in a dictionnary holding the species performance at a given time step.
It runs at initialization, to set up niche partitioning. It used the  with mean temperature for the year `mean_annual` calculated before hand

# Arguments
- `mean_annual`: mean annual temperature for the following year. `x` in the Gaussian function takes its value.
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

"""
    updatefitness!(sppref, mean_opt::Float64, std_tol::Float64, mean_annual::Float64, max_fitness::Float64)

Calculate the fitness value according to a Gauss function:

         f(x) = a*exp(-((x-b)²)/(2*c²)),

and store the value in a dictionnary holding the species performance at a given time step.
Fitness is updated every begining of the year, with mean temperature for the year `mean_annual` calculated before hand.

# Arguments
- `mean_annual`: mean annual temperature for the following year. `x` in the Gaussian function takes its value.
- `maximal fitness::Float64=1.0`: parameter `a`, the height of the curve's peak.
- `mean_opt::Float64`: parameter `b` is the position of the center of the peak.
- `std_tol::Float64`: parameter `c` is the standard deviation,
"""
function updatefitness!(sppref::SppRef, mean_annual::Float64, max_fitness::Float64, t::Int64, settings::Dict{String, Any})
    for sp in sppref.sp_id

        mean_opt = sppref.temp_opt[sp]
        std_tol = sppref.temp_tol[sp]

        absolute_fitness = max_fitness*exp(-((mean_annual-mean_opt)^2)/(2*(std_tol^2)))
        sppref.fitness[sp] = absolute_fitness
    end
    open(joinpath(settings["outputat"], settings["simID"], "spp_fitness.csv"),"a") do fitnessfile
    	writedlm(fitnessfile, hcat(t, sp, get(sppref.fitness, sp, "NA")) for sp in collect(keys(sppref.fitness)))
    end
end

"""
    simulate!()
Run all functions
"""
function simulate()

    #TODO put init into its own function

    settings = parse_commandline()
    println(keys(settings))

    tdist = nothing
    if settings["disturbtype"] != "none"
        tdist = loaddisturbance(settings)
    end
    println(tdist)

    id_counter = 0
    management_counter = 0

    landpars = read_landpars(settings)
    sppref = read_spinput(settings)
    traitranges = define_traitranges(settings)
    interaction, scen, remaining = implicit_insect(settings)
    mylandscape, landavail = submodels.landscape_init(landpars)
    K = initK(landavail, settings, 1)
    T, mean_annual = updateenv!(1, landpars)
    updatefitness!(sppref, mean_annual, 1.0)
    plants, id_counter = initplants(landavail, sppref, id_counter, settings, K)
    
    # check-points
    println("Land init stored in object of type $(typeof(landpars))")
    println("Sp info stored in object of type $(typeof(sppref))")
    println("Sp ranges stored in object of type $(typeof(traitranges))")
    println("Landscape initialized: type $(typeof(mylandscape))")
    println("Plants initialized: type $(typeof(plants))")
    println("Starting simulation")

    # ORGANIZE OUTPUT FOLDERS
    #########################

    srand(settings["rseed"])

    results_folder = string()

    for rep in 1:settings["nreps"]

        #TODO
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
        open(joinpath(simresults_folder, "simsettings.jl"),"w") do ID
            println(ID, "tdist = ", tdist)
            println(ID, "landpars = ", typeof(landpars), "\ninitial = ", typeof(landpars.initialland), "\ndisturb = ", typeof(landpars.disturbland))
            println(ID, "interaction = ", interaction)
            println(ID, "scen = ", scen)
            println(ID, "remaining = ", remaining)
            println(ID, "commandsettings = ", settings)
        end

        # INITIALIZE FILE TO LOG SIMULATION PROGRESS
        open(joinpath(simresults_folder, "checkpoint.txt"),"w") do sim
            println(sim,string("Simulation: ",settings["simID"],now()))
        end

        # INITIALIZE FILE TO LOG LIFE-HISTORY EVENTS
        open(joinpath(simresults_folder, "eventslog.txt"),"w") do sim
            writedlm(sim, hcat("week", "event", "stage", "age"))
        end

        # INITIALIZE FILE TO LOG METABOLIC RATES
        open(joinpath(settings["outputat"],settings["simID"],"metaboliclog.txt"),"a") do sim

            writedlm(sim, hcat("stage", "age", "rate", "probability", "event"))

        end

        # INITIALIZE FILE TO LOG SEED PRODUCTION
        open(joinpath(simresults_folder, "offspringproduction.csv"),"w") do seedfile
            writedlm(seedfile, hcat(["week" "sp" "stage" "mode"], "abundance"))
        end

	    # INITIALIZE FILE TO LOG SPECIES FITNESS
        open(joinpath(simresults_folder, "spp_fitness.csv"),"w") do fitnessfile
            writedlm(fitnessfile, hcat("week", "sp", "fitness"))
        end

        # MODEL RUN
        ############
        for t in 1:settings["timesteps"]

            # check-point
            open(joinpath(joinpath(simresults_folder, "checkpoint.txt")),"a") do sim
                println(sim, "\nWEEK $t")
		        println(sim, "Species richness: $(length(unique(map(x -> x.sp, plants))))")
            end

	        println("\nWEEK $t")
	        println("Species richness: $(length(unique(map(x -> x.sp, plants))))")

            # APPLY LANDSCAPE DISTURBANCE
            if settings["disturbtype"] in ["frag" "loss"] && t in tdist
                landscape, landavail = disturb!(mylandscape,landavail,plants,t,settings,landpars,tdist)
            end
            updateK!(K, landavail, settings, t, tdist)

	        if rem(t, 52) == 1
                T, mean_annual = updateenv!(t, landpars)
	        else
	            T = updateenv!(t, landpars)
	        end

            updatefitness!(sppref, mean_annual, 1.0, t, settings)

            orgstable(plants,t,settings)

            if ((31 < rem(t,52) < 39) && management_counter < 1)
                management_counter = manage!(plants, t, management_counter, settings)
            end

            biomass_production = sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants), 0.00001))
           open(joinpath(simresults_folder, "checkpoint.txt"),"a") do sim
                writedlm(sim, hcat("Biomass production:", biomass_production))
            end

            allocate!(plants, t, aE, Boltz, settings, sppref, T, biomass_production, K, "a")
            survive!(plants, t, cK, K, settings, sppref, landavail, T, biomass_production, "a")
	    allocate!(plants, t, aE, Boltz, settings, sppref, T, biomass_production, K, "j")
            survive!(plants, t, cK, K, settings, sppref, landavail, T, biomass_production, "j")
            develop!(plants, settings, t)
            mate!(plants, t, settings, scen, tdist, remaining)
            id_counter = mkoffspring!(plants, t, settings, sppref, id_counter, landavail, T, traitranges)
            seedsi = getreleases(plants, t)
            disperse!(landavail, seedsi, plants, t, settings, sppref, landpars, tdist)
            establish!(plants, t, settings, sppref, T, biomass_production, K)
            shedd!(plants, sppref, t, settings)
        end
    end
    return settings, results_folder
end

# run simulation
settings, results_folder = simulate()

# analyse results
analysED(settings, results_folder)
