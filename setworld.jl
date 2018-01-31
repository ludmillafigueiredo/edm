module setworld

using Distributions
srand(123) # set seed

export WordCell,

include("temp_prec.jl")
#TODO feed dynamics landscape change

mutable struct WordCell
	suitability::Bool
	resourcesFON_sum::Float64
	pollinationFON_sum::Float64
	temperature::Float64
	precipitation::Float64
	plant::Any
end

# Outer constructors: default values
WordCell() = WordCell(0,0,0,0,0,nothing)

# Functions

function landscape_init(fragmentation::Bool, xlength::Int, ylength::Int, n_frags::Int) #TODO GIS file and control how fragments become unsuitable
	# '''
	# Set the initial world conditions for the number of fragments (n_frags) and the individuals. Landscape is a an array of 3 dimensions: xlength and ylength are maximal dimensions of fragments, n_frags are the number of fragments, when in a fragmentation scenario
	# '''
	if fragmentation == false #TODO see how fragmentatin affects the construction of fragments
		landscape = fragment_init(fragmentation,xlength,ylength)
	else
	# landscape = Array{WordCell}(xlength,ylength, n_frags)
	# fill!(landscape,WordCell())

		for f in 1:n_frags
			landscape[:,:,f] = fragment_init(fragmentation,xlength,ylength)
		end
	end
    return landscape
end

# Unity test: consider
# landtest = landscape_init(false, 5,5, 1)
# # Test
# open("Edlandscape.txt", "a+") do io
# 	dump(landtest)
#  	write(io, (landtest))
# end

function fragment_init!(#landscape::Array{WordCell},
	fragmentation::Bool,
	xlength::Int64,
	ylength::Int64) #TODO read x and y length from external file
	# '''
	# Set environmental conditions on cells of fragments. The fragments are either random or x and y map to a spatial configuration. Initial environmental conditions are read from temp_prec file
	# '''
	# Create basic fragment environmental conditions
	fragment = Array{WordCell}(xlength,ylength)
	fill!(fragment,WordCell())

	# # Create the 4 matrices that keep track of each fragment: #TODO if metacommunity structure is implicit, might not be necessary
	# #TODO should they be a 3d structure like landscape or a sparse matriz? Does global make sense in each of these cases?
	# global res_FONs_matrix = Matrix{Dict}(size(fragment))
	# global poll_FONs_matrix = Matrix{Dict}(size(fragment))
	# global pollinators_matrix = Matrix{Pollinator}(size(fragment)) #TODO it should include herbivores, later

	for y in 1:ylength, x in 1:xlength
		fragment[x,y].suitability = true
		fragment[x,y].temperature = rand(Normal(meantemp,tempsd),1)[1]
		fragment[x,y].precipitation = rand(Normal(meanprec,precsd),1)[1]
		#TODO map different means for different patchs?
	end
	return fragment
end

landscape_init(false, 5,5, 1)

function populate(xlength,ylength, agent::Plant)
	# '''
	# Populate the fragment: randomly draw some xs and ys for placing new individuals. Used for all agents (one method for each).
	# '''
	rng = MersenneTwister(1234) #TODO check this RNG stuff ~ like seed?
	X <- zeros(xlength)
	Y <- zeros(ylength)

	rand!(rng, X, collect(1:xlength))
	rand!(rng, Y, collect(1:ylength))

	for y in Y, x in X
		if fragment[x,y].resourcesFON_sum < 1
			### Plants are stored in fragments matrixes
			fragment[x,y].Plant = createPlant(x,y,res_FONs_matrix) #TODO make sure two plants can't occupy the same location
	end
	return fragment
end

function populate(xlength,ylength, agent::Plant)
### Pollinators are stored in the pollinator matrix
pollinators_matrix[x,y]] = createPollinator(x,y)
