# Functions for reading inputs and initialisation of variables of the model

"""
using RCall
"""

function parse_commandline()
    sets = ArgParseSettings() #object that will be populated with the arguments by the macro
    @add_arg_table! sets begin
        "--simID"
        help = "Name of the folder where outputs will be stored."
        arg_type = String
        default = "test_100_100_100"

        "--rseed"
        help = "Seed for RNG"
        arg_type = Int
        default = 777

        "--sppinput"
        help = "Name of file with species list."
        arg_type = String
		default = "examples/perform_optim/100spp_sppinput.csv"

        "--pollination"
        help = "How to explicitly model insects:
                pollination-independent reproduction \"indep\";
                equal pollination loss for all species \"equal\"."
        arg_type = String
		default = "examples/test/pollination_example.jl"

        "--initial_land"
        help = "Name of file with landscape size values: areas of fragments, mean (and s.d.) temperature."
        arg_type = String
		default = "examples/perform_optim/control_100m2.grd"

        "--disturb_type"
        help = "Type of disturbance to be implemented: none, area_loss, poll_loss, area+poll_loss, clim+area_loss, clim+poll_loss, clim+area+poll_loss"
        arg_type = String
        default = "none"

        "--disturb_land"
        help = "When simulating area loss, the path to a table with disturbance times (td) and proportions of loss of habitat area (template: habitatloss_template.csv). When simulating fragmentation, the path to a table with disturbance times (td) and proportions of loss of habitat area (template: habitatfrag_template.csv)"
        arg_type = Any
        default = nothing

        "--output_freq"
        help = "Frequency of output (number of weeks)"
        arg_type = Int
        default = 26

        "--temp_ts"
        help = "Name of file with weekly temperature  and precipitation time series"
        arg_type = String
		default = "examples/perform_optim/temp_160years.csv"

        "--outputat"
        help = "Name of directory where output should be written ."
        arg_type = String
		default = "."


    end

    return parse_args(sets)

end

"""
read_landpars
"""
function read_landpars(settings)

	if settings["disturb_type"] in ["area_loss", "area+poll_loss"]
                disturbance = CSV.read(settings["disturb_land"], DataFrame, header = true,
			      	       types = Dict("td" => Int64, "proportion" => Float64))
        elseif settings["disturb_type"] == "frag"
                disturbance = CSV.read(settings["disturb_land"], DataFrame, header = true,
			      	       types = Dict("td" => Int64, "frag_file" => String))
        else
		disturbance = nothing
	end

	land_pars = LandPars(settings["initial_land"], disturbance)

    return land_pars

end

"""
    read_sppinput(settings)
Reads in species trait values (min. and max. values for the evolvable ones). Stores them in `sppref`, a structure with parameters as Dictionnary fields, where species names are the keys to the parameter values.
"""
function read_sppinput(settings::Dict{String,Any})

    sppinputtbl = CSV.read(settings["sppinput"], DataFrame)

	if !in("species", names(sppinputtbl))
		speciesfile = split(settings["sppinput"],"_")[1]*"_"*split(settings["sppinput"],"_")[2]*"ids.csv"
		sppinputtspecies = CSV.read(speciesfile, DataFrame)
		sppinputtbl = innerjoin(sppinputtbl, sppinputtspecies, on = :sp_id)
	end

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

    sppinputtbl = CSV.read(settings["sppinput"], DataFrame)

	if !in("species", names(sppinputtbl))
		speciesfile = split(settings["sppinput"],"_")[1]*"_"*split(settings["sppinput"],"_")[2]*"ids.csv"
		sppinputtspecies = CSV.read(speciesfile, DataFrame)
		sppinputtbl = innerjoin(sppinputtbl, sppinputtspecies, on = :sp_id)
	end

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

    if occursin("poll", settings["disturb_type"])
       include(abspath(settings["pollination"]))
       disturbed_regime= CSV.read(pollination_file, DataFrame, header = true,
				  types = Dict("td" => Int64, "remaining" => Float64))
       poll_pars = PollPars(pollination_scen, disturbed_regime)
    else
       poll_pars = PollPars("indep", nothing)
    end

    return poll_pars

end

"""
    set_tdist()
Set disturbance time(s).
"""
function set_tdist(settings)

    # select file according to keyword: loss, frag, temp
    if occursin("none", settings["disturb_type"])
        tdist = nothing
    elseif occursin("area", settings["disturb_type"])
        tdist = CSV.read(settings["disturb_land"], DataFrame, header=true, types=Dict("td"=>Int64))[:td]
    end

    return tdist

end
