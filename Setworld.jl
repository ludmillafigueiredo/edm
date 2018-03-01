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

export WorldCell, landscape_init, neighborhood

#Types
mutable struct WorldCell
	avail::Bool
	#resourcesFONs::Dict #separate Dict matrix
	#pollinationFONs::Dict  separate Dict matrix
	temp::Float64
	precpt::Float64
	neighs::Dict
	WorldCell() = new()
end

# Functions
function landscape_init(fragmentation::Bool, xlength::Int, ylength::Int, n_frags::Int, meantemp::Float64, tempsd::Float64, meanprec::Float64,precsd::Float64) #TODO GIS file and control how fragments become unsuitable
	#TODO if fragmentation is false, no need for n_frags
	# landscape = Array{WorldCell}(xlength,ylength, n_frags)
	# fill!(landscape,WorldCell())
	landscape = []
	#TODO built an outeer constructor with values not filled by the for loop

	fragmentation = false # for now, always false: TODO implement fragmentation

	if fragmentation == false
		for frag in 1:n_frags
			#fragment_init!(landscape[:,:,f],fragmentation,xlength,ylength)
			for y in 1:ylength, x in 1:xlength
				newcell = WorldCell()
				newcell.avail = true
				newcell.temp = rand(Normal(meantemp,tempsd),1)[1]
				newcell.precpt = rand(Normal(meanprec,precsd),1)[1]
				newcell.neighs = Dict()
				push!(landscape,newcell)
			end
		end
		landscape = reshape(landscape, xlength, ylength) #! x inner loops matches reshape order
	else
		# TODO represent fragmentation: will depend on writing a funciton to read txt input: store them in a mutable struct makes it easier to call, wherever size they have
	end
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

#TODO create matrix with Connectivity between fragments
#TODO function to update landscape when changing envoronmental conditions

#TODO
"""
Update landscape with individualsavailability locations.
"""
function updatelandscape!()
end

# Unity test #TODO
# check
# landscape_init(false, 5,5, 1)
# landscape_init(false,5,5,2)
end
