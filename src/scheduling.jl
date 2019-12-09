"""
    simulate!()
Run all functions
"""
function run_scheduling(settings, id_counter, management_counter, landpars, poll_pars, K, T, mean_annual, plants,landscape)

    Random.seed!(settings["rseed"])
    
    for rep in 1:settings["nreps"]

        global simresults_folder, results_folder = orga_outputs()
	log_settings()
	output_sppref(SPP_REF)
	
        # RUN MODEL
        ############
        for t in 1:length(temp_ts.week)

            # check-point
            open(joinpath(joinpath(simresults_folder, "checkpoint.txt")),"a") do sim
                println(sim, "\nWEEK $t")
		println(sim, "Species richness: $(length(unique(getfield.(plants, :sp))))")
            end

	    println("\nWEEK $t")
	    println("Species richness: $(length(unique(getfield.(plants, :sp))))")

            # APPLY LANDSCAPE DISTURBANCE
	    if landpars.disturbance != nothing
                if settings["disturb_type"] in ["frag" "loss"] && t in landpars.disturbance.td
                    landscape = disturb!(landscape,plants,t,settings,landpars)
                end
	    end
	    
            updateK!(K, landscape, settings, t, landpars)

	    if rem(t, 52) == 1
                T, mean_annual = setenv!(t, temp_ts)
	    else
	        T = setenv!(t, temp_ts)
	    end
	    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	        println(sim, "Temperature for week $t: $T")
	        println(sim, "Mean for the year of week $t: $mean_annual")
	    end
	    
            updatefitness!(mean_annual, 1.0, t, settings)

            write_output(plants,t,settings)

            if ((31 < rem(t,52) < 39) && management_counter < 1)
                management_counter = manage!(plants, t, management_counter, settings)
            end

            biomass_production = sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]),plants),NOT_0))
            open(joinpath(simresults_folder, "checkpoint.txt"),"a") do sim
                println(sim, "Biomass production: $biomass_production")
            end

	    die_seeds!(plants, settings, t, T)

	    # Adults growth and mortality
            grow!(plants, t, settings, T, biomass_production, K, "a")
            die!(plants, settings, T, "a", t)
	    compete_die!(plants, t, settings, landscape, T, "a")
	    
	    # Juvenile growth and mortality 
	    grow!(plants, t, settings, T, biomass_production, K, "j")
            die!(plants, settings, T, "j", t)
	    compete_die!(plants, t, settings, landscape, T, "j")

	    # all individuals except seeds in flowers get older over time
	    s_flower = filter(x -> x.stage == "s-in-flower", plants)
	    filter!(x -> x.stage != "s-in-flower", plants)
	    map(x -> age!(x), plants) # surviving get older
	    append!(plants, s_flower)
	    
	    develop!(plants, settings, t)

	    mate!(plants, t, settings, poll_pars::PollPars)

	    # Offspring production
	    id_counter = mkseeds!(plants, settings, id_counter, T, t)
	    id_counter = clone!(plants, settings, id_counter, t)

	    id_counter = self_pollinate!(plants, settings, id_counter, t)
	    
	    setfield!.(plants, :mated, false) # plant only produces seeds again if it gets pollinated
	    
	    justdispersed = disperse!(landscape, plants, t, settings,  landpars)

	    establish!(justdispersed, plants, t, settings,  T, biomass_production, K)
	    
	    shedflower!(plants, t, settings)

	    if (rem(t,52) == 51)
	       winter_dieback!(plants, t)
	    end

	    println("id_counter: $id_counter, # plants: $(length(plants))")
        end
    end

    return settings, results_folder

end
