"""
    This module contains the type of the cell and functions for setting up initial environmental conditions and changing it when necessary.

        WorldCell type
            This type has fields designed to store all relevant environmental conditions and plants steming points. Other types of organisms are sotred in arrays and have a location, but do not need to "point" to a landscape cell.

                landscape_init()
                landscape_init(fragmentation, xlength,ylength, n_frags, meantemp, tempsd, meanprec, precsd)

                landscape_init() creates a multidimensional array composed of habitat grid cells (Array{WorldCell}(xlength, ylength, n_frags)). If -fragmentation` is `false`, then one single grid of `xlength` x `ylength` cells is created. If `fragmentation` is `true`, `n_frags`
                Arguments
                `fragmentation::Bool` statement about fragmentation of the landscape
                `xlength::Int64` Maximal x length of a fragment or total x length of an unfragmented landscape
                `ylength::Int64` Maximal y length of a fragment or total y length of an unfragmented landscape
                `n_frags::Int64` Number of fragments
                `meantemp::Float64` and `tempsd::Float64` mean temperature of a patch and standard deviation
                `meanprec::Float64` and `precsd::Float64` mean precipitation  of a patch and standard-deviation
                """
module Setworld

using Distributions
using Organisms

export LandPars, landscape_init, updateenv!, destroyarea!, fragment!

const tK = 273.15 # °C to K converter

#Store simulation parameters referent to the landscape
mutable struct LandPars
    # Original landscape parameters
    npatches::Int # patches of habitat (fragments are after disturbance)
    plength::Array{Int64,1}
    initialarea::Int64
    initialconnect::Any # might not be simulated.
    # Disturbed landscape parameters 
    nfrags::Any # might not be simulated
    flength::Any # might not be simulated
    disturbarea::Any # might not be simulated 
    disturbconnect::Any # might not be simulated
    # general
    bufferarea::Any # might not be simulated
    meantempts::Array{Float64,1} # all fragments get the same temperature
end

#############
# Functions #
#############
"""
landscape_init()
Create the initial landscape structure.
"""
function landscape_init(landpars::LandPars)

    landscape = Array{Dict{String,Float64}} #if created only inside the loop, remains a local variable
    
    for p in collect(1:landpars.npatches)

        patch = fill(Dict{String,Float64}(), landpars.plength[p],landpars.plength[p])
	
	if p == 1
	    landscape = patch #when empty, landscape cant cat with frag
	else
	    landscape = cat(3,landscape, patch)
	end
    end

    landavail = fill(true,size(landscape))
    
    return landscape, landavail
end

"""
updateenv!(landscape,t)
Update temperature and precipitation values according to the weekly input data (weekly means and ).
"""
function updateenv!(t::Int64, landpars::LandPars)

    T = landpars.meantempts[t] + tK
    #unity test
    println("Temperature for week $t: $T")
    
    return T
end

"""
destroyarea!()
Destroy proportion of habitat area according to input file. Destruction is simulated by making affected cells unavailable for germination and killing organisms in them.
"""
<<<<<<< Updated upstream
function destroyarea!(landmode::String, landpars::LandPars, landavail::Array{Bool,N} where N, loss::Float64, settings::Dict{String,Any})
=======
function destroyarea!(landavail::Array{Bool,3}, landpars::LandPars, settings::Dict{String,Any})
>>>>>>> Stashed changes

    if settings["landmode"] == "artif"
        # DESTROY HABITAT
        # index of the cells still availble:
        available = find(x -> x == true, landavail)
        # number of cells to be destroyed:
        loss = landpars.disturbarea/landpars.initialarea # just a proportion, unit doesn't matter
        lostarea = round(Int,loss*length(available), RoundUp) 

        # Unity test
        if lostarea > length(available)
            error("Destroying more than the available area.")
        end

        # go through landscape indexes of the first n cells from the still available
        for cell in available[1:lostarea] 
            landavail[cell] = false # and destroy them
        end

        #unity test
        open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            println(sim, "Number of destroyed cells: $lostarea")
        end
        
    elseif settings["landmode"] == "real"
        # rebuild the landscape according to shape file
        landscape = Array{Dict{String,Float64}}
        
        for frag in collect(1:landpars.nfrags)

            fragment = fill(Dict{String,Float64}(), landpars.flength[frag],landpars.flength[frag])
	    
	    if frag == 1
	        landscape = fragment #when empty, landscape cant cat with frag
	    else
	        landscape = cat(3,landscape, frag)
	    end
            
        end

        landavail = fill(true,size(landscape))
        
    end
end

"""
    fragment!()
This function is only called for simulating the fragmentation of an originally single landscape. 
"""
function fragment!(landscape::Array{Dict{String,Float64},N} where N, settings::Dict{String,Any}, landpars::LandPars, orgs::Array{Organisms.Organism,1})

    # Built fragmented landscape
    for frag in collect(1:landpars.nfrags)

        fragment = fill(Dict{String,Float64}(), landpars.flength[frag],landpars.flength[frag])
	
	if frag == 1
	    newlandscape = fragment #when empty, landscape cant cat with frag
	else
	    newlandscape = cat(3,landscape, frag)
	end
    end

    # Resettle the individuals in the new landscape:
    # list all indexes of old landscape
    landscapeidx = collect(collect.indices(landscape)) # first collect returns a tuple
    # create array to store landscape cells that were not destroyed
    remaincells = fill(Array{Int64}, size(newlandscape)) # the size 
    # sample old cells to fill the new landscape
    sample!(landscapeidx, remaincells; replace = false)
    # kill the organisms that have not remained in the new landscape configuration
    filter!(x -> x.location in remaincells, orgs)
    # update their .location field
    newloc = []
    for o in orgs
        newloc <- [ind2sub(remaincells,
                           find(x->x == collect(orgs[o].location), remaincells))]
        orgs[o].location = Tuple(newloc)
    end
        
    landscape = newlandscape
    landavail = fill(true,size(landscape))
    
    return landscape, landavail
    
end

"""
disconnect!()
Decreases the connectivity of an already fragmented landscape.
"""
function disconnect!()
    
end

end
