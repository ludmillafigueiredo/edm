"""
This module contains functions related to outputting raw results from the model, as well as default derived analysis.
"""

function orga_outputs()

    results_folder = joinpath(settings["outputat"], settings["simID"])
    
    check_overwrite(results_folder)
    mkpath("$(results_folder)")
    
    # INITIALIZE FILE TO LOG SIMULATION PROGRESS
    open(joinpath(results_folder, "checkpoint.txt"),"w") do sim
        println(sim,string("Simulation: ",settings["simID"],now()))
    end

    # INITIALIZE FILE OF MAIN OUTPUT
    header = hcat(["week"],
                  reshape(collect(string.(fieldnames(Plant)[1:19])),1,:),
                  ["leaves" "stem" "root" "repr"],
                  reshape(collect(string.(fieldnames(Plant)[21:end])),1,:))
    open(joinpath(settings["outputat"],settings["simID"],"statevars_ind.txt"), "w") do output
        writedlm(output, header) #reshape(header, 1, length(header)))
    end

    # INITIALIZE FILE TO OUTPUT SEED PRODUCTION
    open(joinpath(results_folder, "offspringproduction.csv"),"w") do seedfile
        writedlm(seedfile, hcat("week", "sp", "stage", "mode", "abundance"))
    end

    # INITIALIZE FILE TO LOG LIFE-HISTORY EVENTS
    open(joinpath(results_folder, "events.csv"),"w") do sim
        writedlm(sim, hcat("week", "event", "stage", "age", "n_events"))
    end

    # INITIALIZE FILE TO LOG SPECIES FITNESS
    open(joinpath(results_folder, "spp_fitness.csv"),"w") do fitnessfile
        writedlm(fitnessfile, hcat("week", "sp", "fitness"))
    end

    # INITIALIZE FILE TO LOG SPECIES FITNESS
    open(joinpath(results_folder, "pollination_log.csv"),"w") do pollfile
        writedlm(pollfile, hcat("week", "vector", "proportion"))
    end

    return results_folder
end

""" 
    outputorgs(plants,t,settings)
Saves a long format table with the organisms field informations.
"""
function write_output(plants::Array{Plant,1}, t::Int64)

    # only info on juveniles and adults is output
    juvs_adlts = findall(x -> x.stage in ["j", "a"], plants)
    
    
    # output plants info
    for o in juvs_adlts
        open(joinpath(results_folder,"statevars_ind.txt"), "a") do output
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

struct OffspringInfo
    sp::String
    stage::String
    repmode::String
    amount::Int64
end

let offspring_out = OffspringInfo[]
    
    global function gather_offspring!(t, sp, stage, repmode, amount)
        
        push!(offspring_out, OffspringInfo(sp, stage, repmode, amount))
        
    end
    
    global function output_offspring!(t)
        # output it
        open(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv"),"a") do f
            for o in eachindex(offspring_out)
                writedlm(f, hcat(t,
                                 offspring_out[o].sp,
                                 offspring_out[o].stage,
                                 offspring_out[o].repmode,
                                 offspring_out[o].amount))
            end
  	end
        # reset it
        offspring_out = OffspringInfo[]
    end
end

"""
analsyED()
Run R script of analysis after the model finishes the simulation
"""
function analysED(settings, land_pars, poll_pars)

    simID = settings["simID"]
    @rput simID
    @rput EDMdir
    disturbance = settings["disturb_type"]
    @rput disturbance
    
    if settings["disturb_type"] in ["area_loss", "area+poll_loss"]
        tdist = land_pars.disturbance.td
    elseif settings["disturb_type"] == "poll"
        tdist = poll_pars.regime.td
    else
	tdist = nothing
    end
    @rput tdist

    outputsdir = joinpath(settings["outputat"],settings["simID"])
    @rput outputsdir

    # run analysis
    R"source(\"src/analysED.R\")"

end
