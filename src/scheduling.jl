"""
    simulate!()
Run all functions
"""
function run_scheduling(settings, tdist, id_counter, management_counter, landpars, interaction, scen, remaining, K, T, mean_annual, plants)

    Random.seed!(settings["rseed"])
    
    for rep in 1:settings["nreps"]

        global simresults_folder, results_folder = orga_outputs()
	log_settings()
	output_sppref(sppref)
	
        # RUN MODEL
        ############
        for t in 1:settings["timesteps"]

            # check-point
            open(joinpath(joinpath(simresults_folder, "checkpoint.txt")),"a") do sim
                println(sim, "\nWEEK $t")
		        println(sim, "Species richness: $(length(unique(getfield.(plants, :sp))))")
            end

	        println("\nWEEK $t")
	        println("Species richness: $(length(unique(getfield.(plants, :sp))))")

            # APPLY LANDSCAPE DISTURBANCE
            if settings["disturbtype"] in ["frag" "loss"] && t in tdist
                global mylandscape, landscape = disturb!(mylandscape,landavail,plants,t,settings,landpars,tdist)
            end
            updateK!(K, landavail, settings, t, tdist)

	        if rem(t, 52) == 1
                T, mean_annual = setenv!(t, landpars)
	        else
	            T = setenv!(t, landpars)
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
            allocate!(plants, t, aE, Boltz, settings, T, biomass_production, K, "a")
            die!(plants, settings, T, "a")
	    compete_die!(plants, t, cK, settings, landavail, T, "a")
	    
	    # Juvenile growth and mortality 
	    allocate!(plants, t, aE, Boltz, settings, T, biomass_production, K, "j")
            die!(plants, settings, T, "j")
	    compete_die!(plants, t, cK, settings, landavail, T, "j")

	    # all individuals except seeds in flowers get older over time
	    s_flower = filter(x -> x.stage == "s-in-flower", plants)
	    filter!(x -> x.stage != "s-in-flower", plants)
	    setfield!.(plants, :age, +(1)) # surviving get older
	    append!(plants, s_flower)
	    
	    develop!(plants, settings, t)

	    mate!(plants, t, settings, scen, tdist, remaining)

	    # Offspring production
	    id_counter = mkseeds!(plants, settings, id_counter, T, t)
	    id_counter = clone!(plants, settings, id_counter)

	    setfield!.(plants, :mated, false) # plant only produces seeds again if it gets pollinated
	    
	    justdispersed = disperse!(landavail, plants, t, settings,  landpars, tdist)

	    establish!(justdispersed, plants, t, settings,  T, biomass_production, K)
	    
	    shedflower!(plants, t, settings)

	    if (rem(t,52) == 51)
	       winter_dieback!(plants, t)
	    end
        end
    end

    return settings, results_folder

end