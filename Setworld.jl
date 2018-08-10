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

export LandPars, landscape_init, updateenv!, destroyarea!

const tK = 273.15 # °C to K converter

#Simulation parameters storage:
mutable struct LandPars
fxlength::Array{Int64,1}
fylength::Array{Int64,1}
meantempts::Array{Float64,1} #all fragments get the same temperature
sdtempts::Array{Float64,1}
meanprects::Array{Float64,1}
sdprects::Array{Float64,1} # this is probably gonna go
nfrags::Int64
end

# Functions
function landscape_init(landpars::LandPars)

    landscape = Array{Dict{String,Float64}} #if created only inside the loop, remains a local variable
    
    for frag in collect(1:landpars.nfrags)

        fragment = fill(Dict{String,Float64}(), landpars.fxlength[frag],landpars.fylength[frag])
	
	if frag == 1
	    landscape = fragment #when empty, landscape cant cat with frag
	else
	    landscape = cat(3,landscape, frag)
	end
    end

    landavail = fill(true,size(landscape))
    
    return landscape, landavail
end

"""
            updateenv!(landscape,t)
            Update temperature and precipitation values according to the weekly input data (weekly means and ).
            """
function updateenv!(landscape::Array{Dict{String,Float64},2}, t::Int64, landpars::LandPars)

    T = landpars.meantempts[t] + tK
    #unity test
    println("Temperature for week $t: $T")

    if t != 1
        fill!(landscape, Dict())
        #occupied = find(x -> !isempty(x), landscape)
        #println("Cells to be reset: $occupied.")
        #for cell in occupied
	    #empty!(landscape[occupied])
        #end
    end
    
    return T
end

"""
            destroyarea!()
            Destroy proportion of habitat area according to input file. Destruction is simulated by making affected cells unavailable for germination and killing organisms in them.
                """
function destroyarea!(landavail::Array{Bool,2}, loss::Float64, settings::Dict{String,Any})
    # DESTROY HABITAT
    # index of the cells still availble:
    available = find(x -> x == true, landavail)
    # number of cells to be destroyed:
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
end

"""
                fragmentation!()
                """
function fragment!(landscape::Array{Dict{String,Float64}}, t, settings::Dict{String,Any})
end

end
