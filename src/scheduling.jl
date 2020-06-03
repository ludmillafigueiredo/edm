"""
    simulate!()
Run all functions
"""
function run_scheduling(settings, management_counter, land_pars, poll_pars, K, T, mean_annual, plants,landscape)

    Random.seed!(settings["rseed"])
    
    
        global results_folder = orga_outputs()
	log_settings()
	output_sppref(SPP_REF)
	
        # RUN MODEL
        ############
        for t in 1:length(temp_ts.week)

            println("\nWEEK $t")
	    println("Species richness: $(length(unique(getfield.(plants, :sp))))")

            # APPLY LANDSCAPE DISTURBANCE
	    if settings["disturb_type"] in ["frag", "area_loss", "area+poll_loss"] &&
	       t in land_pars.disturbance.td
                landscape = disturb!(landscape,plants,t,settings,land_pars)
	    end
	    
            updateK!(K, landscape, settings, t, land_pars)

	    if rem(t, 52) == 1
                T, mean_annual = setenv!(t, temp_ts)
	    else
	        T = setenv!(t, temp_ts)
	    end

            updatefitness!(mean_annual, 1.0, t, settings)

            write_output(plants,t,settings)

            if ((31 < rem(t,52) < 39) && management_counter < 1)
                management_counter = manage!(plants, t, management_counter, settings)
            end

            biomass_production = sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]),plants),NOT_0))

            check_ages(plants)
	    
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
	    
	    mature!(plants, settings, t)

	    check_ages(plants)
	    
	    mate!(plants, t, settings, poll_pars::PollPars)

	    # Offspring production
	    mkseeds!(plants, settings, T, t)
	    clone!(plants, settings, t)

	    self_pollinate!(plants, settings, t)
	    
	    setfield!.(plants, :mated, false) # plant only produces seeds again if it gets pollinated
	    
	    println("Dispersal:")
	    @time disperse!(landscape, plants, t, settings,  land_pars)

	    println("Dispersal:")
	    @time establish!(plants, t, settings,  T, biomass_production, K)
	    
	    shedflower!(plants, t, settings)

	    if (rem(t,52) == 51)
	       winter_dieback!(plants, t)
	    end
	    
        end
end
