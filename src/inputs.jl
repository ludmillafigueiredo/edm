# Functions for reading inputs and initialisation of variables of the model

using RCall

function parse_commandline()
    sets = ArgParseSettings() #object that will be populated with the arguments by the macro
    @add_arg_table sets begin
        "--simID"
        help = "Name of the folder where outputs will be stored."
        arg_type = String
        default = "profiling"

        "--nreps"
        help = "Number of replicates."
        arg_type = Int64
        default = 1

        "--rseed"
        help = "Seed for RNG"
        arg_type = Int
        default = 100

        "--outputat"
        help = "Name of directory where output should be written ."
        arg_type = String
        default = "outputs"

        "--spinput"
        help = "Name of file with species list."
        arg_type = String
        default = joinpath("test_inputs", "weisssingles/weissplants_sppinput.csv")

        "--insect"
        help = "How to explicitly model insects:
                pollination-independent reproduction \"indep\";
                equal pollination loss for all species \"equal\"."
        arg_type = String
        default = joinpath("test_inputs", "insects_indep.csv")

        "--initialland"
        help = "Name of file with landscape size values: areas of fragments, mean (and s.d.) temperature."
        arg_type = String
        default = joinpath("test_inputs","control_49m2.grd")

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
        default = joinpath("inputs", "landbuffer.jl")

        "--disturbland"
        help = "Either a shape file (if \`landmode\`)"
        arg_type = Any
        default = nothing

        "--timesteps"
        help = "Duration of simulation in weeks."
        arg_type = Int
        default = 5200

        "--tout"
        help = "Frequency of output (number of weeks)"
        arg_type = Int
        default = 12

        "--temp_ts"
        help = "Name of file with weekly temperature  and precipitation time series"
        arg_type = String
        default = joinpath("test_inputs", "temp1917_2017.csv")

        "--timemsg"
        help = "Output timing to terminal, as well as checkpoint"
        arg_type = Bool
        default = false
    end
    return parse_args(sets)
end

"""
    read_landpars(settings)
Read, derive and store values related to the environmental conditions: temperature, landscape size, availability, and changes thereof.
"""
function read_landpars(settings::Dict{String,Any})
    # Temperature time-series
    temp_tsinput = CSV.read(settings["temp_ts"])

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
        landconfig = rcopy(Array{Any}, reval(joinpath(EDMdir, "landconfig.R")))
        # store them in landpars
        # in the absence of disturbance-related files, R returns Nullable values
        landpars = LandPars(landconfig[[1]],
                                     areatocell(landconfig[[2]]),
                                     landconfig[[3]],
                                     landconfig[[4]],
                                     landconfig[[5]],
                                     areatocell(landconfig[[6]]),
                                     landconfig[[7]],
                                     landconfig[[8]],
                                     temp_tsinput[:,:meantemp])

    elseif settings["landmode"] == "artif" # simulation arena is built from a neutral landscape

        if settings["disturbtype"] in ["frag" "loss"]

            # send file names to R
            initialland = settings["initialland"]
            @rput initialland

            disturbland = settings["disturbland"] #object has to be initialized
            if settings["disturbtype"] == "loss"
                disturbland = CSV.read(settings["disturbland"], header = true, types = Dict("td" => Int64, "proportion" => Float64))
            elseif settings["disturbtype"] == "frag"
                disturbland = CSV.read(settings["disturbland"], header = true, types = Dict("td" => Int64, "disturbland" => String))
            end
            @rput disturbland

            # get 0-1 habitat suitability matrices
            landconfig = rcopy(Array{Any}, reval(string("source(\"",joinpath(EDMdir,"landconfig.R"),"\")")))
            @rget initialmatrix
            @rget disturbmatrix

            # if fragmentation is simulated, the file names are assgined to dist field, if area loss, the proportion holding columns, if none, default value nothing
            if disturbmatrix == "notfrag"
                landpars = NeutralLandPars(Int.(initialmatrix),
                                                    disturbland,
                                                    temp_tsinput[:,:meantemp])
            else
                landpars = NeutralLandPars(Int.(initialmatrix),
                                                    Int.(disturbmatrix),
                                                    temp_tsinput[:,:meantemp])
            end
            #TODO block above can handle pollination
        else
            # send file names to R
            initialland = settings["initialland"]
            @rput initialland
            disturbland = settings["disturbland"]
            @rput disturbland

            # get suitability matrices from R
            landconfig = rcopy(Array{Any}, reval(string("source(\"",joinpath(EDMdir, "landnlmconfig.R"),"\")")))
            @rget initialmatrix
            @rget disturbmatrix

            landpars = NeutralLandPars(Int.(initialmatrix),
                                                nothing,
                                                temp_tsinput[:,:meantemp])
        end
    else
        error("Please inform how landscape should be simulated:\n
               artificially created (artif)\n
               built from .grd input files (real).")
    end
    return landpars
end

"""
    read_sppinput(settings)
Reads in species trait values (min. and max. values for the evolvable ones). Stores them in `sppref`, a structure with parameters as Dictionnary fields, where species names are the keys to the parameter values.
"""
function read_sppinput(settings::Dict{String,Any})

    spinputtbl = CSV.read(settings["spinput"])

    sppref = SppRef()

    # sp_id, clonality and fitness columns cannot be set as the rest
    sppref.sp_id = spinputtbl[:, :sp_id]
    sppref.clonality = Dict(spinputtbl[:, :sp_id][i] => spinputtbl[:, :clonality][i] == "true"
		        	      		    for i in 1:length(spinputtbl[:, :sp_id]))
    sppref.kernel = Dict(spinputtbl[:, :sp_id][i] => spinputtbl[:, :kernel][i]
		        	      		    for i in 1:length(spinputtbl[:, :sp_id]))
    sppref.fitness = Dict()

    println(fieldnames(typeof(sppref))[4:end-1])
    for field in fieldnames(typeof(sppref))[4:end-1]
    	setfield!(sppref, field,
	          Dict(spinputtbl[:, :sp_id][i] => Float64(spinputtbl[:, field][i])
		        	      		    for i in 1:length(spinputtbl[:, :sp_id])))
    end
    
    return sppref
end

"""
    define_traitranges(settings)
Based on input species traits values, set the limits upon which they can vary during sexual reproduction: [min., max.].
Values are store in an object of type TraitRanges, where each field is a Dictionary for a trait, where the keys are the species, and the values are Arrays holding min. and max. values the trait can take.
"""
function define_traitranges(settings::Dict{String,Any})

    spinputtbl = CSV.read(settings["spinput"])

    traitranges = TraitRanges(
        Dict(spinputtbl[:,:sp_id][i] =>
             [spinputtbl[:,:compartsize][i], spinputtbl[:,:compartsize][i]]
             for i in 1:length(spinputtbl[:,:sp_id])),
        Dict(spinputtbl[:,:sp_id][i] =>
             [Int(round(spinputtbl[:,:span_min][i], RoundDown)), Int(round(spinputtbl[:,:span_max][i], RoundUp))]
             for i in 1:length(spinputtbl[:,:sp_id])),
        Dict(spinputtbl[:,:sp_id][i] =>
             [Int(round(spinputtbl[:,:firstflower_min][i], RoundDown)), Int(round(spinputtbl[:,:firstflower_max][i], RoundUp))]
             for i in 1:length(spinputtbl[:,:sp_id])),
        Dict(spinputtbl[:,:sp_id][i] =>
             [Int(round(spinputtbl[:,:floron][i], RoundDown)), Int(round(spinputtbl[:,:floron][i], RoundUp))]
             for i in 1:length(spinputtbl[:,:sp_id])),
        Dict(spinputtbl[:,:sp_id][i] =>
             [Int(round(spinputtbl[:,:floroff][i], RoundDown)), Int(round(spinputtbl[:,:floroff][i], RoundUp))]
             for i in 1:length(spinputtbl[:,:sp_id])),
        Dict(spinputtbl[:,:sp_id][i] =>
             [Int(round(spinputtbl[:,:seednumber_min][i], RoundDown)), Int(round(spinputtbl[:,:seednumber_max][i], RoundUp))]
             for i in 1:length(spinputtbl[:,:sp_id])),
        Dict(spinputtbl[:,:sp_id][i] =>
             [Int(round(spinputtbl[:,:seedon][i], RoundDown)), Int(round(spinputtbl[:,:seedon][i], RoundUp))]
             for i in 1:length(spinputtbl[:,:sp_id])),
        Dict(spinputtbl[:,:sp_id][i] =>
             [Int(round(spinputtbl[:,:seedoff][i], RoundDown)), Int(round(spinputtbl[:,:seedoff][i], RoundUp))]
             for i in 1:length(spinputtbl[:,:sp_id])),
        Dict(spinputtbl[:,:sp_id][i] =>
             [Int(round(spinputtbl[:,:bankduration_min][i], RoundDown)), Int(round(spinputtbl[:,:bankduration_max][i], RoundUp))]
             for i in 1:length(spinputtbl[:,:sp_id])))
    return traitranges
end

"""
    read_insects(settings)
Reads how insects are going to be implicitly simulated.
"""
function read_insects(settings::Dict{String,Any})

    insectsinput = CSV.read(settings["insect"])
    interaction = insectsinput[:, :interaction][1] # type of interaction which is affected
    scen = insectsinput[:, :scen][1] # pollination scenario
    remaining = insectsinput[:, :remaining] # proportion of pollination service loss
    return interaction, scen, remaining
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
            tdist = CSV.read(settings["insect"])[:, :td]
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

