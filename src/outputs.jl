"""
This module contains functions related to outputting raw results from the model, as well as default derived analysis.
"""

# Julia packages
using RCall
using DelimitedFiles

function orga_outputs()

        results_folder = joinpath(settings["outputat"], settings["simID"])

	try
            mkpath("$(results_folder)")
            println("Output will be written to '$(results_folder)'")
        catch
            println("Overwriting results to existing '$(results_folder)' folder")
        end

        # INITIALIZE FILE TO LOG SIMULATION PROGRESS
        open(joinpath(results_folder, "checkpoint.txt"),"w") do sim
            println(sim,string("Simulation: ",settings["simID"],now()))
        end

        # INITIALIZE FILE TO LOG LIFE-HISTORY EVENTS
        open(joinpath(results_folder, "eventslog.txt"),"w") do sim
            writedlm(sim, hcat("week", "event", "stage", "age"))
        end

        # INITIALIZE FILE TO LOG METABOLIC RATES
        open(joinpath(settings["outputat"],settings["simID"],"metaboliclog.txt"),"w") do sim
            writedlm(sim, hcat("stage", "age", "rate", "probability", "event"))
        end

        # INITIALIZE FILE TO LOG SEED PRODUCTION
        open(joinpath(results_folder, "offspringproduction.csv"),"w") do seedfile
            writedlm(seedfile, hcat(["week" "sp" "stage" "mode"], "abundance"))
        end

	    # INITIALIZE FILE TO LOG SPECIES FITNESS
        open(joinpath(results_folder, "spp_fitness.csv"),"w") do fitnessfile
            writedlm(fitnessfile, hcat("week", "sp", "fitness"))
        end

	return results_folder
end

function output_sppref(SPP_REF)
	 	# printout SPP_REF
	open(joinpath(results_folder, "sppref_traitvalues.csv"), "w") do ref
            writedlm(ref, reshape(collect(string.(fieldnames(SppRef))), 1,:), ",")
	end
    
	for sp in SPP_REF.sp
    	    sp_ref = reshape(collect(map(x -> getfield(SPP_REF,x)[sp],fieldnames(SppRef)[2:end])),1,:)
    	    open(joinpath(results_folder, "sppref_traitvalues.csv"), "a") do ref
                writedlm(ref, [sp sp_ref], ",")
            end
	end
end

function log_settings()
        open(joinpath(results_folder, "simsettings.jl"),"w") do ID
            println(ID, "landpars = $(repr(typeof(landpars))) \ninitial = $(repr(typeof(landpars.initial))) \ndisturb_land = $(repr(typeof(landpars.disturbance))) \ndisturb_poll = $(repr(poll_pars))")
            println(ID, "commandsettings = $(repr(settings))")
        end
end

"""
    outputorgs(plants,t,settings)
Saves a long format table with the organisms field informations.
"""
function write_output(plants::Array{Plant,1}, t::Int64, settings::Dict{String,Any})

    # only info on juveniles and adults is output
    juvs_adlts = findall(x -> x.stage in ["j", "a"], plants)
    
    # output header
    if t == 1
        header = hcat(["week"],
                      reshape(collect(string.(fieldnames(Plant)[1:19])),1,:),
                      ["leaves" "stem" "root" "repr"],
                      reshape(collect(string.(fieldnames(Plant)[21:end])),1,:))
        open(joinpath(settings["outputat"],settings["simID"],"statevars_ind.txt"), "w") do output
            writedlm(output, header) #reshape(header, 1, length(header)))
        end
    end

    # output plants info
    if t == 1 || rem(t,settings["tout"]) == 0

        for o in juvs_adlts
            open(joinpath(settings["outputat"],settings["simID"],"statevars_ind.txt"), "a") do output
                writedlm(output, hcat(t,
                                      plants[o].id,
                                      plants[o].stage,
                                      plants[o].location,
                                      plants[o].sp,
                                      plants[o].kernel,
                                      plants[o].clonality,
				      plants[o].pollen_vector,
				      plants[o].self_failoutcross,
				      plants[o].self_proba,
                                      plants[o].compartsize,
                                      plants[o].span,
                                      plants[o].firstflower,
                                      plants[o].floron,
                                      plants[o].floroff,
                                      plants[o].seednumber,
                                      plants[o].seedon,
                                      plants[o].seedoff,
                                      plants[o].bankduration,
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
log_abovegroundmass(when, proc)
"""
function log_abovegroundmass(when,proc)

    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        println(sim,"Above-grd mass $when $proc: $(sum(vcat(map(x -> sum(values(x.mass))-x.mass["root"],
	                                                        filter(x -> x.stage in ["j", "a"], plants)), NOT_0)))g")
    end
    
end

function log_sim(msg)

    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        println(sim, msg)
    end
    
end

function log_age()

    open(joinpath(results_folder, "checkpoint.txt"),"a") do sim
        println(sim, "Juvs age >= firstflower: $(filter(x -> (x.stage == "j" && x.age > x.firstflower), 
                                                  plants) |> x -> [getfield.(x, :id), 
                                                                   getfield.(x, :stage), 
                                                                   getfield.(x, :firstflower), 
                                                                   getfield.(x, :age)])")
    end
    
end
"""
analsyED()
Run R script of analysis after the model finishes the simulation
"""
function analysED(settings)

    parentsimID = settings["simID"]
    @rput parentsimID
    EDdir = pwd()
    @rput EDdir
    nreps = settings["nreps"]
    @rput nreps
    disturbance = settings["disturb_type"]
    @rput disturbance

    outputsdir = results_folder
    @rput outputsdir

    # run analysis
    R"source(\"src/analysED.R\")"

end