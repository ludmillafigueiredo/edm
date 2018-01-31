module setworld
# Contains functions for

using Distributions
srand(123) # set seed

export WordCell, meantemp, tempsd, meanprec,

include("temp_prec.jl") #TODO parse?
#TODO feed dynamics landscape change

mutable struct WordCell
	suitability::Bool
	resourcesFONs::Dict
	pollinationFONs::Dict
	resourcesFON_sum::Float64
	pollinationFON_sum::Float64
	temperature::Float64
	precipitation::Float64
	plant::Any
end

# Outer constructors: default values
WordCell() = WordCell(0,Dict(),Dict(),0,0,0,0,nothing)

# Functions

function landscape_init(fragmentation::Bool, xlength::Int, ylength::Int, n_frags::Int) #TODO GIS file and control how fragments become unsuitable
	# '''
	# Set the initial world conditions for the number of fragments (n_frags) and the individuals. Landscape is a an array of 3 dimensions: xlength and ylength are maximal dimensions of fragments, n_frags are the number of fragments, when in a fragmentation scenario
	# '''
	# if fragmentation == false #TODO see how fragmentatin affects the construction of fragments
	# 	landscape = fragment_init(fragmentation,xlength,ylength)
	# else
	# create basic empty landscape
	landscape = Array{WordCell}(xlength,ylength, n_frags)
	fill!(landscape,WordCell())

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
		# TODO represent fragmentation
	end
    return landscape
end

#TODO create matrix with Connectivity between fragments

# Unity test #TODO
# landtest = landscape_init(false, 5,5, 1)
# # Test
# open("Edlandscape.txt", "a+") do io
# 	dump(landtest)
#  	write(io, (landtest))
# end
landscape_init(false, 5,5, 1)
landscape_init(false,5,5,2)
