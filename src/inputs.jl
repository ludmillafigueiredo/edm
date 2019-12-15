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

        "--sppinput"
        help = "Name of file with species list."
        arg_type = String
        default = joinpath("test_inputs", "weisssingles/weissplants_sppinput.csv")

        "--pollination"
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
        help = "Type of environmental disturbance to be implemented: habitat area loss \"area_loss\", pollination loss \"poll_loss\", area and pollination loss \"area+poll_loss\" or temperature change \"temp\""
        arg_type = String
        default = "none"

        "--disturb_file"
        help = "When simulating area loss, the path to a table with disturbance times (td) and proportions of loss of habitat area (template: habitatloss_template.csv). When simulating fragmentation, the path to a table with disturbance times (td) and proportions of loss of habitat area (template: habitatfrag_template.csv)"
        arg_type = Any
        default = nothing

        "--tout"
        help = "Frequency of output (number of weeks)"
        arg_type = Int
        default = 12

        "--temp_ts"
        help = "Name of file with weekly temperature  and precipitation time series"
        arg_type = String
        default = joinpath("template_files", "temperaturegoettingen_18572017.csv")

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
	    
	if settings["disturb_type"] in ["area_loss", "area+poll_loss"]
                disturbance = CSV.read(settings["disturb_file"], header = true,
			      	       types = Dict("td" => Int64, "proportion" => Float64))
        elseif settings["disturb_type"] == "frag"
                disturbance = CSV.read(settings["disturb_file"], header = true,
			      	       types = Dict("td" => Int64, "frag_file" => String))	    
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

    sppinputtbl = CSV.read(settings["sppinput"])

    sppref = SppRef()

    # species, clonality and fitness columns cannot be set as the rest
    sppref.sp = sppinputtbl[:, :species]
    sppref.clonality = Dict(sppinputtbl[:, :species][i] => sppinputtbl[:, :clonality][i] == "TRUE"
		        	      		    for i in 1:length(sppinputtbl[:, :species]))
    sppref.kernel = Dict(sppinputtbl[:, :species][i] => sppinputtbl[:, :kernel][i]
		        	      		    for i in 1:length(sppinputtbl[:, :species]))
    sppref.pollen_vector = Dict(sppinputtbl[:, :species][i] => sppinputtbl[:, :pollen_vector][i]
		        	      		    for i in 1:length(sppinputtbl[:, :species]))
    sppref.self_failoutcross = Dict(sppinputtbl[:,:species][i] => sppinputtbl[:,:self_failoutcross][i] == "TRUE"
		        	      		    for i in 1:length(sppinputtbl[:, :species]))
    sppref.fitness = Dict()

    for field in fieldnames(typeof(sppref))[6:end-1]
    	setfield!(sppref, field,
	          Dict(sppinputtbl[:, :species][i] => Float64(sppinputtbl[:, field][i])
		        	      		    for i in 1:length(sppinputtbl[:, :species])))
    end
    
    return sppref
end

"""
    define_traitranges(settings)
Based on input species traits values, set the limits upon which they can vary during sexual reproduction: [min., max.].
Values are store in an object of type TraitRanges, where each field is a Dictionary for a trait, where the keys are the species, and the values are Arrays holding min. and max. values the trait can take.
"""
function define_traitranges(settings::Dict{String,Any})

    sppinputtbl = CSV.read(settings["sppinput"])

    traitranges = TraitRanges(
        Dict(sppinputtbl[:,:species][i] =>
             [sppinputtbl[:,:compartsize][i], sppinputtbl[:,:compartsize][i]]
             for i in 1:length(sppinputtbl[:,:species])),
        Dict(sppinputtbl[:,:species][i] =>
             [Int(round(sppinputtbl[:,:span_min][i], RoundDown)), Int(round(sppinputtbl[:,:span_max][i], RoundUp))]
             for i in 1:length(sppinputtbl[:,:species])),
        Dict(sppinputtbl[:,:species][i] =>
             [Int(round(sppinputtbl[:,:firstflower_min][i], RoundDown)), Int(round(sppinputtbl[:,:firstflower_max][i], RoundUp))]
             for i in 1:length(sppinputtbl[:,:species])),
        Dict(sppinputtbl[:,:species][i] =>
             [Int(round(sppinputtbl[:,:floron][i], RoundDown)), Int(round(sppinputtbl[:,:floron][i], RoundUp))]
             for i in 1:length(sppinputtbl[:,:species])),
        Dict(sppinputtbl[:,:species][i] =>
             [Int(round(sppinputtbl[:,:floroff][i], RoundDown)), Int(round(sppinputtbl[:,:floroff][i], RoundUp))]
             for i in 1:length(sppinputtbl[:,:species])),
        Dict(sppinputtbl[:,:species][i] =>
             [Int(round(sppinputtbl[:,:seednumber_min][i], RoundDown)), Int(round(sppinputtbl[:,:seednumber_max][i], RoundUp))]
             for i in 1:length(sppinputtbl[:,:species])),
        Dict(sppinputtbl[:,:species][i] =>
             [Int(round(sppinputtbl[:,:seedon][i], RoundDown)), Int(round(sppinputtbl[:,:seedon][i], RoundUp))]
             for i in 1:length(sppinputtbl[:,:species])),
        Dict(sppinputtbl[:,:species][i] =>
             [Int(round(sppinputtbl[:,:seedoff][i], RoundDown)), Int(round(sppinputtbl[:,:seedoff][i], RoundUp))]
             for i in 1:length(sppinputtbl[:,:species])),
        Dict(sppinputtbl[:,:species][i] =>
             [Int(round(sppinputtbl[:,:bankduration_min][i], RoundDown)), Int(round(sppinputtbl[:,:bankduration_max][i], RoundUp))]
             for i in 1:length(sppinputtbl[:,:species])))
    return traitranges
end

"""
    read_pollination(settings)
Reads how insects are going to be implicitly simulated.
"""
function read_pollination(settings::Dict{String,Any})

    if settings["disturb_type"] in ["poll_loss", "area+poll_loss"]
        poll_pars = PollPars("indep", nothing)
    else
       include(settings["pollination"])
       disturbed_regime= CSV.read(pollination_file, header = true,
				  types = Dict("td" => Int64, "remaining" => Float64))
       poll_pars = PollPars(pollination_scen, disturbed_regime)	
    end
    
    return poll_pars

end

"""
    set_tdist()
Set disturbance time(s).
"""
function set_tdist(settings)

    # select file according to keyword: loss, frag, temp
        if settings["disturb_type"] == "none"
            tdist = nothing
	elseif settings["disturb_type"] in ["area_loss" "area+poll_loss"]
            tdist = CSV.read(settings["disturb_file"], header=true, types=Dict("td"=>Int64))[:td]
        elseif settings["disturb_type"] == "temp"
            println("Temperature change is simulated with the temperature file provided.")
        end

	return tdist
end
