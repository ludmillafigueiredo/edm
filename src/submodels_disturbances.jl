"""
    disturb!(landscape, landscape, plants, t, settings, landpars, tdist)
Change landscape structure and metrics according to the simulation scenario (`loss`, habitat loss, or `frag`, fragmentation)
"""
function disturb!(landscape::BitArray{2}, plants::Array{Plant,1}, t::Int64, settings::Dict{String,Any}, landpars::LandPars)

    if settings["disturb_type"] in ["area_loss", "area+poll_loss"]
         landscape = destroyarea!(landpars, landscape, settings, t)
    elseif settings["disturb_type"] == "frag"
         landscape = load_fragmentation!(landpars, landscape, settings, t)
    end

    destroyorgs!(plants, landscape, settings)

    return landscape
end

"""
destroyarea!()
Destroy proportion of habitat area according to input file. Destruction is simulated by making affected cells unavailable for germination and killing organisms in them.
"""
function destroyarea!(landpars::LandPars, landscape::BitArray{2}, settings::Dict{String,Any}, t::Int64)

    # index of the cells still available, 
    available = findall(x -> x == true, landscape)
    # number of cells to be destroyed:
    loss = landpars.disturbance[landpars.disturbance.td .== t, :proportion][1]	

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
function load_fragmentation!(landpars::LandPars, landscape::BitArray{2}, settings::Dict{String,Any}, t::Int64)

     disturb_file = landpars.disturbance[landpars.disturbance.td .== t, :disturb_file][1]
     @rput disturb_file
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

"""
disturb_pollinate!()
Call the pollinate!() function for groups of individuals according to 
"""
function disturb_pollinate!(flowering::Array{Plant,1}, poll_pars::PollPars, plants::Array{Plant,1}, t::Int64)

    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Pollination disturbed by $(poll_pars.scen) ...")
    end

    insects_plants = filter(x -> occursin("insects", x.pollen_vector), flowering)
    total_p = 0
    
    if poll_pars.scen  == "equal" # all plants species are equally affected
        for sp in unique(getfield.(insects_plants, :sp))
	    npoll = pollinate!(insects_plants, sp, flowering, plants, poll_pars, t)
	    total_p += npoll
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

    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Percentage of pollinated flowers under disturbance $(total_p/length(insects_plants)) ...")
    end

    
end