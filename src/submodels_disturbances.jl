"""
    disturb!(landscape, landscape, plants, t, settings, land_pars, tdist)
Change landscape structure and metrics according to the simulation scenario (`loss`, habitat loss, or `frag`, fragmentation)
"""
function disturb!(landscape::BitArray{2}, plants::Array{Plant,1}, t::Int64, settings::Settings, land_pars::LandPars)

    if settings.disturb_type in ["area_loss", "area+poll_loss"]
         landscape = destroyarea!(land_pars, landscape, settings, t)
    elseif settings.disturb_type == "frag"
         landscape = load_fragmentation!(land_pars, landscape, t)
    end

    destroyorgs!(plants, landscape, settings)

    return landscape
end

"""
destroyarea!()
Destroy proportion of habitat area according to input file. Destruction is simulated by making affected cells unavailable for germination and killing organisms in them.
"""
function destroyarea!(land_pars::LandPars, landscape::BitArray{2}, settings::Settings, t::Int64)

    # index of the cells still available,
    available = findall(x -> x == true, landscape)
    # number of cells to be destroyed:
    loss = land_pars.disturbance[land_pars.disturbance.td .== t, :proportion][1]

    lostarea = round(Int,loss*length(available), RoundUp)

    # unit test
    if lostarea > length(available)
        error("Destroying more than the available area.")
    end

    # go through landscape indexes of the first n cells from the still available
    for cell in available[1:lostarea]
        landscape[cell] = false # and destroy them
    end

    open(abspath(joinpath(settings.outputat,settings.simID,"checkpoint.txt")),"a") do sim
        println(sim, "Remaining available area: $(findall(x -> x == true, landscape)) m2")
    end

    return landscape

end


"""
fragment!()
Create new landscape from raster file of a frafmented one
"""
function load_fragmentation!(land_pars::LandPars, landscape::BitArray{2}, t::Int64)

     disturb_land = land_pars.disturbance[land_pars.disturbance.td .== t, :disturb_land][1]
     @rput disturb_land
     R"source(\"src/load_fragmentation.R\")"
     @rget disturb_matrix

     # convert raster matrix to BitArray (smaller than Bool)
     landscape = BitArray(disturb_matrix)

     return landscape

end

"""
    destroyorgs!(plants)
Kill organisms that where in the lost habitat cells.
"""
function destroyorgs!(plants::Array{Plant,1}, landscape::BitArray{2}, settings::Settings)

    kills = []
    for o in 1:length(plants)
	if landscape[plants[o].location...] == false
	    push!(kills,o)
	end
    end

    if length(kills) > 0 # trying to delete at index 0 generates an error
	deleteat!(plants, kills)
    end

    # check-point
    open(abspath(joinpath(settings.outputat,settings.simID,"checkpoint.txt")),"a") do sim
        println(sim, "Killed $(length(kills))")
    end

end

# method for disturb_poll
function get_npoll(insects_plants::Array{Plant,1}, poll_pars::PollPars, t::Int64)

    npoll_dflt = rand(Distributions.Binomial(Int(ceil(length(insects_plants)*VST_DFLT)),
						 INSCT_EFFC))[1]
	if poll_pars.scen == "indep"
     	    npoll = npoll_dflt
	elseif poll_pars.scen in ["equal", "rdm"]
	    if t in poll_pars.regime.td
	        npoll = Int(ceil(npoll_dflt * poll_pars.regime[poll_pars.regime.td.== t,
		                                               :remaining][1]))
	    else
		npoll = npoll_dflt
	    end
	elseif poll_pars.scen == "spec"
	    # pseudo
	end

    return npoll

end

"""
disturb_pollinate!()
Call the pollinate!() function for groups of individuals according to
"""
function disturb_pollinate!(flowering::Array{Plant,1}, poll_pars::PollPars, plants::Array{Plant,1}, t::Int64, settings::Settings)

    # check-point
    open(joinpath(settings.outputat,settings.simID,"checkpoint.txt"),"a") do sim
        println(sim, "Pollination disturbed by $(poll_pars.scen) ...")
    end

    insects_plants = filter(x -> occursin("insects", x.pollen_vector), flowering)
    total_p = 0

    if poll_pars.scen  == "equal" # all plants species are equally affected
        for sp in unique(getfield.(insects_plants, :sp))
	    insects_plants_sp = filter(x -> x.sp == sp, insects_plants)
    	    npoll = get_npoll(insects_plants_sp, poll_pars, t)
	    log_pollination(flowering, npoll, "insects-dist", t)
            pollinated = sample(insects_plants_sp, npoll, replace = false,
                                ordered = false)
            filter!(x -> !(x.id in getfield.(pollinated, :id)), flowering)
            setfield!.(pollinated, :mated, true)
            append!(plants, pollinated)
            total_p += npoll
            # unit test
            check_duplicates(plants)
	end
    elseif poll_pars.scen == "rdm"
    	npoll = get_npoll(insects_plants, poll_pars, t)
	total_p += npoll
	if npoll > 0
	    pollinate!(insects_plants, plants, npoll)
	end
    elseif poll_pars.scen == "spec"

    else
	error("Please chose a pollination scenario \"scen\" in insect.csv:
                       - \"indep\": sexual reproduction happens independently of pollination
                       - \"rmd\": random loss of pollinator species
                       - \"spec\": specific loss of pollinator species")

    end

end
