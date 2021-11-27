function log_initialabund(plants)

    plants_spp = getfield(plants, :sp)
    initial_abunds = Dict(i => sum(plants_spp .== i) for i in unique(plants_spp))

    open(joinpath(results_folder, "initialabundances.txt"),"w") do sim
        println(sim, initial_abunds)
    end

end

function log_sppref(SPP_REF)
    # prinoutput_freq SPP_REF
    open(joinpath(results_folder, "sppref_traitvalues.csv"), "w") do ref
        writedlm(ref, reshape(collect(string.(fieldnames(SppRef))), 1,:), ",")
        for sp in SPP_REF.sp
    	    sp_ref = reshape(collect(map(x -> getfield(SPP_REF,x)[sp],fieldnames(SppRef)[2:end])),1,:)
    	    writedlm(ref, [sp sp_ref], ",")
        end
    end
end

function log_settings()
        open(joinpath(results_folder, "simsettings.jl"),"w") do ID
            println(ID, "land_pars = $(repr(typeof(land_pars))) \ninitial = $(repr(typeof(land_pars.initial))) \ndisturb_land = $(repr(typeof(land_pars.disturbance))) \ndisturb_poll = $(repr(poll_pars))")
            println(ID, "commandsettings = $(repr(settings))")
        end
end

struct EventInfo
    event::String
    stage::String
    mean_age::Float64
    amount::Int64
end

let events_log = EventInfo[]

    global function gather_event!(t, event, stage, mean_age, amount)

        push!(events_log, EventInfo(event, stage, mean_age, amount))

    end

    global function log_events!(t)
        # output it
        open(joinpath(results_folder,"events.csv"),"a") do f
            for e in eachindex(events_log)
                writedlm(f, hcat(t,
                                 events_log[e].event,
                                 events_log[e].stage,
                                 events_log[e].mean_age,
                                 events_log[e].amount))
            end
  	end
        # reset it
        events_log = EventInfo[]
    end
end

"""
log_abovegroundmass(when, proc)
"""
function log_abovegroundmass(when,proc)

    open(abspath(joinpath(settings.outputat, settings.simID,"checkpoint.txt")),"a") do sim
        println(sim,"Above-grd mass $when $proc: $(sum(vcat(map(x -> sum(values(x.mass))-x.mass["root"],
	                                                        filter(x -> x.stage in ["j", "a"], plants)), NOT_0)))g")
    end

end

function log_sim(msg)

    open(abspath(joinpath(settings.outputat, settings.simID,"checkpoint.txt")),"a") do sim
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

function log_pollination(flowering, npoll, pollen_vector, t)

    open(joinpath(results_folder, "pollination_log.csv"),"a") do pollfile
        writedlm(pollfile, hcat(t, pollen_vector, npoll/(length(flowering)+NOT_0)))
    end

end
