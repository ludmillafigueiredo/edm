"""
    simulate!()
Run all functions
"""
function run_scheduling(settings, management_counter, land_pars, poll_pars, K, T, mean_annual, plants,landscape)

    Random.seed!(settings.rseed)


    global results_folder = orga_outputs()
    log_settings()
    log_sppref(SPP_REF)

    # RUN MODEL
    ############
    for t in 1:length(temp_ts.week)

        println("\nWEEK $t")
	println("Species richness: $(length(unique(getfield.(plants, :sp))))")

        # APPLY LANDSCAPE DISTURBANCE
	if settings.disturb_type in ["frag", "area_loss", "area+poll_loss"] &&
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

        if t == 1 || rem(t,settings.output_freq) == 0
            write_output(plants,t)
            output_offspring!(t)
            log_events!(t)
        end

        if ((31 < rem(t,52) < 39) && management_counter < 1)
            management_counter = manage!(plants, t, management_counter)
        end

        biomass_production = sum(vcat(map(x -> (x.mass.leaves+x.mass.stem),plants),NOT_0))

	print("check_ages:")
    @time check_ages(plants)

	print("die_seeds:")
	@time die_seeds!(plants, settings, t, T)

	# Adults growth and mortality
	print("A grow:")
    @time grow!(plants, t, T, biomass_production, K, "a")
	print("A die:")
    @time die!(plants, settings, T, "a", t)
	print("A compete_die:")
	@time compete_die!(plants, t, landscape, T, "a")

	# Juvenile growth and mortality
	print("J grow:")
	@time grow!(plants, t, T, biomass_production, K, "j")
	print("J die:")
    @time die!(plants, settings, T, "j", t)
	print("J compete_die:")
	@time compete_die!(plants, t, landscape, T, "j")

	# all individuals except seeds in flowers get older over time
	s_flower = filter(x -> x.stage == "s-in-flower", plants)
	filter!(x -> x.stage != "s-in-flower", plants)
	map(x -> age!(x), plants) # surviving get older
	append!(plants, s_flower)

	print("mature:")
	@time mature!(plants, t)

	print("check_ages:")
	@time check_ages(plants)

	print("mate:")
	@time mate!(plants, t, poll_pars::PollPars, settings)

	# Offspring production
	mkseeds!(plants, settings, T, t)
	clone!(plants, settings, t)

	self_pollinate!(plants, settings, t)

	setfield!.(plants, :mated, false) # plant only produces seeds again if it gets pollinated

	print("Dispersal:")
	@time disperse!(landscape, plants, t, settings,  land_pars)

	print("Establish:")
	@time establish!(plants, t, settings,  T, biomass_production, K)

	shedflower!(plants, t)

	if (rem(t,52) == 51)
	    winter_dieback!(plants, t)
	end

    end
end
