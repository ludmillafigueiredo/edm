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

        "--initial_land"
        help = "Name of file with landscape size values: areas of fragments, mean (and s.d.) temperature."
        arg_type = String
        default = joinpath("test_inputs","control_49m2.grd")

        "--disturb_type"
        help = "Type of environmental disturbance to be implemented: habitat area loss \"loss\", habitat fragmentation \"frag\" or temperature change \"temp\""
        arg_type = String
        default = "none"

        "--disturb_file"
        help = "When simulating area loss, the path to a table with disturbance times (td) and proportions of loss of habitat area (template: habitatloss_template.csv). When simulating fragmentation, the path to a table with disturbance times (td) and proportions of loss of habitat area (template: habitatfrag_template.csv)"
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
read_landpars
"""
function read_landpars(settings)
	    
	if settings["disturb_type"] == "loss"
                disturbance = CSV.read(settings["disturb_file"], header = true,
			      	       types = Dict("td" => Int64, "proportion" => Float64))
        elseif settings["disturb_type"] == "frag"
                disturbance = CSV.read(settings["disturb_file"], header = true,
			      	       types = Dict("td" => Int64, "disturb_file" => String))	    
        else
		disturbance = nothing
	end

	landpars = LandPars(settings["initial_land"], disturbance)

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
    set_tdist()
Set disturbance time(s).
"""
function set_tdist(settings)

    # select file according to keyword: loss, frag, temp
        if settings["disturb_type"] == "none"
            tdist = nothing
	elseif settings["disturb_type"] in ["loss" "frag"]
            tdist = CSV.read(settings["disturb_file"], header=true, types=Dict("td"=>Int64))[:td]
        elseif settings["disturb_type"] == "poll"
            tdist = CSV.read(settings["insect"])[:, :td]
        elseif settings["disturb_type"] == "temp"
            println("Temperature change is simulated with the temperature file provided.")
        else
            error("Please specify one of the disturbance scenarios with `--disturb`:
                \n\'none\' if no disturbance should be simulated,
                \n\'loss\' for habitat area loss,
                \n\'frag\' for habitat fragmentation,
                \n\'temp\' for temperature change,
                \n\'poll\' for pollination loss.")
        end

	return tdist
end
