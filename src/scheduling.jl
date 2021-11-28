"""
    simulate!()
Run all functions
"""
function run_scheduling(settings, management_counter, land_pars, poll_pars, K, T, mean_annual, landscape)

    Random.seed!(settings.rseed)


    global results_folder = orga_outputs()
    log_settings()
    log_sppref(SPP_REF)

    # RUN MODEL
    ############
    for t in 1:length(temp_ts.week)

        println("\nWEEK $t")

	#build a list of all plants in the landscape at the start of the timestep, mostly used for debug messages.
	#THIS LIST DOES NOT UPDATE DURING THE TIMESTEP!!!
	plantlist = Plant[]
	for plant_vector in landscape.plants
		plantlist = vcat(plantlist, plant_vector)
	end

	println("Species richness: $(length(unique(getfield.(plantlist, :sp))))")
	println("Total population: ", length(plantlist))

        # APPLY LANDSCAPE DISTURBANCE
	if settings.disturb_type in ["frag", "area_loss", "area+poll_loss"] &&
	    t in land_pars.disturbance.td
            landscape = disturb!(landscape,t,settings,land_pars)
	end

        updateK!(K, landscape, settings, t, land_pars)

	if rem(t, 52) == 1
            T, mean_annual = setenv!(t, temp_ts)
	else
	    T = setenv!(t, temp_ts)
	end

        updatefitness!(mean_annual, 1.0, t, settings)
        if t == 1 || rem(t,settings.output_freq) == 0
            write_output(plantlist,t)
            output_offspring!(t)
            log_events!(t)
        end

        if ((31 < rem(t,52) < 39) && management_counter < 1)
            management_counter = manage_matrix!(landscape.plants, t, management_counter)
        end

        biomass_production = sum(vcat(map(x -> (x.mass.leaves+x.mass.stem),plantlist),NOT_0))

	print("check_ages:")
    @time check_ages_matrix!(landscape.plants)

	print("die_seeds:")
	@time die_seeds_matrix!(landscape.plants, settings, t, T)

	# Adults growth and mortality
	print("A grow:")
    @time grow_matrix!(landscape.plants, t, T, biomass_production, K, "a")
	print("A die:")
    @time die_matrix!(landscape.plants, settings, T, "a", t)
	print("A compete_die:")
	@time compete_die_matrix!(landscape.plants, t, landscape.habitability, T, "a")

	# Juvenile growth and mortality
	print("J grow:")
	@time grow_matrix!(landscape.plants, t, T, biomass_production, K, "j")
	print("J die:")
    @time die_matrix!(landscape.plants, settings, T, "j", t)
	print("J compete_die:")
	@time compete_die_matrix!(landscape.plants, t, landscape.habitability, T, "j")

	print("age_plants:")
	@time age_plants_matrix!(landscape.plants)

	print("mature:")
	@time mature_matrix!(landscape.plants, t)

	print("check_ages:")
	@time check_ages_matrix!(landscape.plants)

	print("mate:")
	@time mate_matrix!(landscape.plants, t, poll_pars, settings)

	# Offspring production
	print("mkseeds:")
	@time mkseeds_matrix!(landscape.plants, settings, T, t)
	print("clone:")
	@time clone_matrix!(landscape.plants, settings, t)

	print("self_pollinate:")
	@time self_pollinate_matrix!(landscape.plants, settings, t)

	# plant only produces seeds again if it gets pollinated
	reset_plant_mating_matrix!(landscape.plants)

	# dispersal is a two step process, first the dispersal matrix is calculated, second the dispersal matrix is applied to the plant matrix
	print("Dispersal Matrix Calculation:")
	@time calculate_dispersal_matrix!(landscape, t, settings,  land_pars)
	#Apply and reset the dispersal matrix
	print("Dispersal Matrix Application:")
	@time apply_dispersal_matrix!(landscape)

	print("Establish:")
	@time establish_matrix!(landscape.plants, t, settings,  T, biomass_production, K)

	shedflower_matrix!(landscape.plants, t)

	if (rem(t,52) == 51)
	    winter_dieback_matrix!(landscape.plants, t)
	end

    end
end
