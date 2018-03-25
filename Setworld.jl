landscape_init
"""
This module contains the type of the cell and functions for setting up initial environmental conditions and changing it when necessary.

# WorldCell type
This type has fields designed to store all relevant environmental conditions and plants steming points. Other types of organisms are sotred in arrays and have a location, but do not need to "point" to a landscape cell.

# landscape_init()
	landscape_init(fragmentation, xlength,ylength, n_frags, meantemp, tempsd, meanprec, precsd)

landscape_init() creates a multidimensional array composed of habitat grid cells (Array{WorldCell}(xlength, ylength, n_frags)). If `fragmentation` is `false`, then one single grid of `xlength` x `ylength` cells is created. If `fragmentation` is `true`, `n_frags`
# Arguments
- `fragmentation::Bool` statement about fragmentation of the landscape
`xlength::Int64` Maximal x length of a fragment or total x length of an unfragmented landscape
`ylength::Int64` Maximal y length of a fragment or total y length of an unfragmented landscape
`n_frags::Int64` Number of fragments
`meantemp::Float64` and `tempsd::Float64` mean temperature of a patch and standard deviation
`meanprec::Float64` and `precsd::Float64` mean precipitation  of a patch and standard-deviation
"""
module Setworld

using Distributions

export WorldCell, landscape_init

#Types
mutable struct WorldCell
	avail::Bool
	#resourcesFONs::Dict #separate Dict matrix
	#pollinationFONs::Dict  separate Dict matrix
	temp::Float64
	precpt::Float64
	neighs::Dict{String,Int64} #
	#stem::Bool #steming point of plant
	#WorldCell() = new() #if this is included, and WorldCell object must be initialized as WordCell() and then completed
end

mutable struct PollCell
	floralres::Dict() #ind => amount of floral resource projected in the cell
end

#WorldCell() = WorldCell(false, 0.0, 0.0, Dict())

# Functions
function landscape_init(simpars)

	landscape = WorldCell[]

	for frag in 1:simpars.nfrags
		for y in 1:simpars.fylength[frag], x in 1:simpars.fxlength[frag]
			newcell = WorldCell(true,
								rand(Normal(simpars.fmeantemp[frag],simpars.ftempsd[frag]),1)[1],
								rand(Normal(simpars.fmeanprec[frag],simpars.fprecsd[frag]),1)[1],
								Dict(),
								false)
			push!(landscape,newcell)
		end
	end

	landscape = reshape(landscape, (simpars.fxlength[1], simpars.fylength[1], simpars.nfrags)) # reshaping is easier then going through every index of a 3D landscap, creating a WordCell cell there and parameterizing it. x in inner loops matches reshape order
	#TODO check reshaping for landscape with frags of different sizes

	return landscape
end

# function read_envcond()
# 	# creat incomplete calls (outer constructors section of Julia manual for when read_envcond(false), conditions don't change)
# end

# """
# """
# function neighborhood(landscape::Array{WorldCell, N} where N)
# 	neighborhood = []
# 	fill!(neighborhood, Dict())
# 	reshape(neighboorhood)
# 	return neighborhood
# end

#TODO create/read in matrix with Connectivity between fragments

"""
Update landscape with individuals availability locations.
"""
function updatelandscape!()
end

#TODO Unity test

#module end!
end
