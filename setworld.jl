"""
This module contains the type of the cell and functions for setting up initial environmental conditions and changing it when necessary.

# WorldCell type
This type has fields designed to store all relevant environmental conditions and plants steming points. Other types of organisms are sotred in arrays and have a location, but do not need to "point" to a landscape cell.

# landscape_init()
Sets up the initial landscape (single patch or fragments) and the initial conditions
Ex:
```julia-repl
```
"""
module setworld

using Distributions

export WordCell, landscape_init

#Types
mutable struct WorldCell
	suitability::Bool
	#resourcesFONs::Dict #separate Dict matrix
	#pollinationFONs::Dict  separate Dict matrix
	temperature::Float64
	precipitation::Float64
	plant::Bool
	resourcesFON_sum::Float64 # orient herbivores
	pollinationFON_sum::Float64 # orient pollinators
	WorldCell() = new()
end

# Functions
function landscape_init(fragmentation::Bool, xlength::Int, ylength::Int, n_frags::Int, meantemp::Float64, tempsd::Float64, meanprec::Float64,precsd::Float64) #TODO GIS file and control how fragments become unsuitable
	landscape = Array{WorldCell}(xlength,ylength, n_frags)
	fill!(landscape,WorldCell())

	fragmentation = false # for now, always TODO still dont know how o implement fragmentation

	if fragmentation == false
		for frag in 1:n_frags
			#fragment_init!(landscape[:,:,f],fragmentation,xlength,ylength)
			for y in 1:ylength, x in 1:xlength
				landscape[x,y,frag].suitability = true
				landscape[x,y,frag].temperature = rand(Normal(meantemp,tempsd),1)[1] #TODO whole fragment is getting the same values
				landscape[x,y,frag].precipitation = rand(Normal(meanprec,precsd),1)[1]
				#TODO map different means for different patchs (change the input file)
			end
		end
	else
		# TODO represent fragmentation: will depend on writing a funciton to read txt input
	end
    return landscape
end

# function read_envcond()
# 	# creat incomplete calls (outer constructors section of Julia manual for when read_envcond(false), conditions don't change)
# end

#TODO create matrix with Connectivity between fragments
#TODO function to change environmental conditions


# Unity test #TODO
# check
# landscape_init(false, 5,5, 1)
# landscape_init(false,5,5,2)
end
