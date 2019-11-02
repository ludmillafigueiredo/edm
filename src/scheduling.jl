"""
    simulate!()
Run all functions
"""
function run_scheduling(settings, tdist, id_counter, management_counter, landpars, sppref, traitranges, interaction, scen, remaining, K, T, mean_annual, plants)

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

            updatefitness!(sppref, mean_annual, 1.0, t, settings)

            write_output(plants,t,settings)

            if ((31 < rem(t,52) < 39) && management_counter < 1)
                management_counter = manage!(plants, t, management_counter, settings)
            end

            biomass_production = sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants), 0.00001))
            open(joinpath(simresults_folder, "checkpoint.txt"),"a") do sim
                writedlm(sim, hcat("Biomass production:", biomass_production))
            end

	    survive!(plants, t, settings, sppref, T)
	    
            allocate!(plants, t, aE, Boltz, settings, sppref, T, biomass_production, K, "a")
            survive!(plants, t, cK, K, settings, sppref, landavail, T, biomass_production, "a")

	    allocate!(plants, t, aE, Boltz, settings, sppref, T, biomass_production, K, "j")
            survive!(plants, t, cK, K, settings, sppref, landavail, T, biomass_production, "j")

	    develop!(plants, settings, t)

	    mate!(plants, t, settings, scen, tdist, remaining, sppref)

	    id_counter = mkoffspring!(plants, t, settings, sppref, id_counter, landavail, T, traitranges)

	    justdispersed = disperse!(landavail, plants, t, settings, sppref, landpars, tdist)

	    establish!(justdispersed, plants, t, settings, sppref, T, biomass_production, K)

	    shedflower!(plants, sppref, t, settings)

	    if (rem(t,52) == 51)
	       winter_dieback!(plants, t)
	    end
        end
    end

    return settings, results_folder

end