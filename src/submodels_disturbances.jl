"""
    disturb!(landscape, landscape, plants, t, settings, land_pars, tdist)
Change landscape structure and metrics according to the simulation scenario (`loss`, habitat loss, or `frag`, fragmentation)
"""
function disturb!(landscape::BitArray{2}, plants::Array{Plant,1}, t::Int64, settings::Dict{String,Any}, land_pars::LandPars)

    if settings["disturb_type"] in ["area_loss", "area+poll_loss"]
         landscape = destroyarea!(land_pars, landscape, settings, t)
    elseif settings["disturb_type"] == "frag"
         landscape = load_fragmentation!(land_pars, landscape, settings, t)
    end

    destroyorgs!(plants, landscape, settings)

    return landscape
end

"""
destroyarea!()
Destroy proportion of habitat area according to input file. Destruction is simulated by making affected cells unavailable for germination and killing organisms in them.
"""
function destroyarea!(land_pars::LandPars, landscape::BitArray{2}, settings::Dict{String,Any}, t::Int64)

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

    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        println(sim, "Remaining available area: $(findall(x -> x == true, landscape)) m2")
    end
	
    return landscape
    
end


"""
fragment!()
Create new landscape from raster file of a frafmented one
"""
function load_fragmentation!(land_pars::LandPars, landscape::BitArray{2}, settings::Dict{String,Any}, t::Int64)

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
function destroyorgs!(plants::Array{Plant,1}, landscape::BitArray{2}, settings::Dict{String,Any})

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
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        println(sim, "Killed $(length(kills))")
    end
    
end

# method for disturb_poll
function get_npoll(tobepoll::Array{Plant,1}, poll_pars::PollPars, t::Int64)

    npoll_dflt = rand(Distributions.Binomial(Int(ceil(length(tobepoll)*VST_DFLT)),
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
pollinate!()
Method used when pollination has been disturbed
"""
function pollinate!(tobepoll_idxs::Array{Plant,1}, npoll::Int64, flowering::Array{Plant,1}, plants::Array{Plant,1}, poll_pars::PollPars, t::Int64)

    log_pollination(npoll, "insects-dist", t)

    poll_idxs = sample(tobepoll_idxs, npoll, replace = false,
                        ordered = false)
    for i in poll_idxs
        plants[i].mated = true
    end
    
    # unit test
    check_duplicates(plants)

end

"""
    disturb_pollinate!()
Call the pollinate!() function for groups of individuals according to 
"""
function disturb_pollinate!(flwr_ids::Array{String,1}, poll_pars::PollPars, plants::Array{Plant,1}, t::Int64)

    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Pollination disturbed by $(poll_pars.scen) ...")
    end

    insct_plnts = filter(x -> occursin("insects", x.pollen_vector) && x.id in flwr_ids, plants)
    
    if poll_pars.scen  == "equal" # allspecies are equally affected
        for sp in unique(getfield.(insct_plants, :sp))
            tobepoll_idxs = findall(x -> x.sp == sp, insct_plants)
	    npoll = get_npoll(tobepoll_idxs, poll_pars, t)
	    pollinate!(tobepoll_idxs, npoll, plants, poll_pars, t)
	end
    elseif poll_pars.scen == "rdm"
        insctplants_idxs = findall(x -> x.id in getfiled.(insct_plants, :id), plants)
        npoll = get_npoll(insctplants_idxs, poll_pars, t)
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
