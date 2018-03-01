#!/usr/bin/env julia
#ID

# chmod u+x <file> make it executable with ./file
# also check chmod744 <file>: 3 type of permission level: user, group and other - 744 refers to each one of them; also true for ls

__precompile()__ #__precompile creates precompile version to reduce loading time
#TODO check the cache stuff for dependencies

module EDmodel

# Load packages
using DataFrames # for reading data
# using Distributions
# using DifferentialEquations

export #types & functions
Plant, Pollinator, WordCell, Seed,

# import?
# The import keyword supports all the same syntax as using, but only operates on a single name at a time. It does not add modules to be searched the way using does. import also differs from using in that functions must be imported using import to be extended with new methods. To do this, import the base of the fct and rewrite it new one with new method(Base.namefunction). Use importall for all functions in a model.

# ---------------------------------------------
#				Composition of types
# ---------------------------------------------
mutable struct Plant
	id::String
	sp::String
	PFT::String
	genotype::Any #TODO own type? traits depend on them
	#TODO traits: traits depend on genotype
	vegmass::Float64
	floral_mass::Float64
	poll_offer::String #temp: should match the pollinator's preference in "resouce_preference"
	pollen_disp_synd::String #wind dispersal or pollinator
	disp_par::Tuple #kernel parameters
	location::Tuple # fixed, cant change
	c_FON::Float64
	#TODO generation size
end

mutable struct Pollinator
	id::String
	sp::String
	genotype::String
	#TODO traits depend on genotype
	biomass::Float64
	disp_par::Float64
	#TODO should the pollinator object hold an interaction matrix or verify it when polinizing? (the latter, but how)
	poll_access::String #TODO
	location::Array{Int,2} # changes at each time step
end

mutable struct WordCell
	suitability::Bool #TODO only describes if this is a calcareous grassland or not. Main reason is the future possibility of feeding GIS that would be read as "grasslands" - suitable - or not. Possibly not relevant for experiment of matrix permeability
	resourcesFON_sum::Float64
	pollinationFON_sum::Float64
	#limiting_availability::Float64 # amount of limiting nutrient available (for plants) #TODO include it
	temperature::Float64
	precipitation::Float64
	#plants_id_list::Set #Which individuals compete in that area easier to check for which other plant individuals have FONs in that cell #TODO should the world landscape hold this info or is the matriz of competition neighborhoods enough?
	plant::Any # or type Any, because it must hold something if not occupied
end

mutable struct Seed
	id::String
	sp::String
	size::Float64 #TODO mass?
end

# Outer constructors for default values
Plant() = Plant(nothing,nothing,nothing,nothing,0,0,nothing,nothing,(0,0),(0,0),0)
Pollinator() = Pollinator(nothing,nothing,nothing,0,0,nothing,[0,0])
WordCell() = WordCell(0,0,0,0,0,nothing)

# Constants
const Boltz = 8.617e-5# Boltzmann constant eV/K (non-SI) 1.38064852e-23 J/K if SI
const Ea = 0.69 # activation energy kJ/mol (non-SI), 0.63eV (MTE - Brown et al. 2004)
const plant_growthrate = exp(25.2) # plant biomass production (Ernest et al. 2003)
#TODO reproduction
#TODO allocation to floral growth: if metabolic rate is the rate of allocation to maintenance, reproduction and growth, then, how would it be modelled? Does allocation to

# ---------------------------------------------
#				    Functions
# ---------------------------------------------
function landscape_init(fragments_file)
	#TODO make it possible to choose fragmetns size
	'''
	Set the initial world conditions for the number of fragments (n_frags) and the individuals.
	'''
	#TODO read fragments file
		#TODO xlength, ylength are the biggest values
		#TODO n_frags is the number of fragments
	#Create landscape 3D array to store fragments #TODO or sparse matrix?
	landscape = Array{WordCell}(xlength,ylength,n_frags) #TODO since fragments can have different sizes, the landscape should be able to store the biggest

	for f in 1:n_frags
		landscape[:,:,f] = fragment_init(fragxlength,fragxlength)
	end
    return landscape
end

function fragment_init(xlength, ylength) #TODO read x and y length from external file
	'''
	Set individiuals and inital environmental conditions on a fragment.
	'''
	# Create basic fragment environmental conditions
	fragment = Array{WordCell}(xlength,ylength)
	fill!(fragment,WordCell())

	# Create the 4 matrices that keep track of each fragment: #TODO if metacommunity structure is implicit, might not be necessary
	#TODO should they be a 3d structure like landscape or a sparse matriz? Does global make sense in each of these cases?
	global res_FONs_matrix = Matrix{Dict}(size(fragment))
	global poll_FONs_matrix = Matrix{Dict}(size(fragment))
	global pollinators_matrix = Matrix{Pollinator}(size(fragment)) #TODO it should include herbivores, later

	for y in 1:ylength, x in 1:xlength
		fragment[x,y].temperature = 25 #TODO read from env data: lines are fragments, columns are niche dimensions (xlength, ylength, temperature,precipitation, umidity, limiting nutrients, etc)
		fragment[x,y].precipitation = 500 #TODO read from env data
	end

	#TODO Populate the fragment: randomly draw some xs and ys for placing new individuals
	rng = MersenneTwister(1234) #TODO check this RNG stuff
	X <- zeros(30)
	Y <- zeros(30)
	#TODO create initial individuals with some resemblance to real abundance or wait for equilibrium?
	rand!(rng, X, collect(1:xlength))
	rand!(rng, Y, collect(1:ylength))

	for y in Y, x in X
		#if fragment[x,y].resourcesFON_sum < 1 #TODO not necessary for initializing but might be for reproduction

		### Plants are stored in fragments matrixes
		fragment[x,y].Plant = createPlant(x,y,res_FONs_matrix) #TODO make sure two plants can't occupay the same location

		### Pollinators are stored in the pollinator matrix
		pollinators_matrix[x,y]] = createPollinator(x,y)
	end
	return fragment
end

function fragment_init!(landscape::Array{WordCell},
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
		fragment[x,y].temperature = rand(Normal(meantemp,tempsd),1)[1] #TODO whole fragment is getting the same values
		fragment[x,y].precipitation = rand(Normal(meanprec,precsd),1)[1]
		#TODO map different means for different patchs?
	end
	return fragment
end

function createPlant(x::Int64,y::Int64)
	# 1. Create sp acording to species data: #TODO use different methods to generate from seeds
	#TODO Read species file
	global plant_id_counter = Int64(1) #TODO after being initialized, it will have to keep track of all other new individuals

	Plant(
	id = string(plant_id_counter),
	sp = "sp1", #TODO read from file
	PFT = "ptf1",#TODO read from file
	genotype = "AA", #TODO own type? traits depend on them
	#TODO traits: traits depend on genotype
	vegmass = 0.1, #TODO read from file
	floral_mass = 0, #TODO read from file
	poll_offer = 'flower',#TODO read from file
	pollen_disp_synd = 'poll' ,#TODO read from file: wind dispersal or pollinator
	disp_par = [1 1], #TODO read from file, kernel parameters
	location = (x,y), # fixed, cant change
	c_FON = 1#TODO make it dependent on biomass
	)
	# Project the neighborhood
	resourceFONproj(x,y,res_FONs_matrix, Plant)

	return Plant
end

function resourceFONproj(x,y,res_FONs_matrix, Plant)
	#TODO there's got to be a way for better moore neighborhood
	#TODO check that cell already has FON sum = 1, no establishment there
	push!(res_FONs_matrix[x,y],Plant.id => 1) #TODO if this multiple atribution will work
	push!(res_FONs_matrix[[x-1,x,x,x+1],[y,y+1,y,y-1]], Plant.id => exp(-fragment[x,y].Plant.c_FON))
	push!(res_FONs_matrix[[x],y], Plant.id => exp(-(fragment[x,y].Plant.c_FON * sqrt(2)))
	#TODO deal with boundaries!!!
end

function createPollinator(x,y)
	global pollinator_id_counter =  Int64(1)

	Pollinator(
	id = string(pollinator_id_counter),
	sp = 'sp1',
	genotype = 'AA',
	#TODO traits depend on genotype
	biomass = 0.01,
	disp_par = (1,1),
	#TODO should the pollinator object hold an interaction matrix or verify it when polinizing? (the latter, but how)
	poll_access = 'flower',
	location = (x,y)
	)
end

function plant_growth(Plant)
	grown_mass += growthrate * Plant.vegmass^(3/4) * exp(-Ea/(Boltz*fragment[Plant.location].temperature)) #TODO same growthrate for insects and plants?
	#TODO can this reference to the fragment be a problem? how will it know which fragment (only  a problem if using the 3D array arrangement for landscape, actually)

	Plant.vegmass += grown_mass

	#TODO growth of floral mass
	#TODO read on minimum biomass for floral allocation
	if plant.biomass > #TODO minimum
		#TODO allocate
end

function pollinator_growth(pollinator::Pollinator)
	grown_mass += growthrate * pollinator.biomass^(3/4) * exp(-Ea/(Boltz*fragment[pollinator.location].temperature)) #TODO same growthrate for insects and plants?
	#TODO can this reference to the fragment be a problem? how will it know which fragment (only  a problem if using the 3D array arrangement for landscape, actually)

	pollinator.biomass += grown_mass
end

function floralFONproj()
end

function pollinator_dispersal() #TODO make an interaction function and use differnt methods depending on partner?
	#1. Check if there is pollination FON where pollinator is
end

function pollination()
end

function plant_reproduction()
	# give
end

function pollinator_reproduction()
end

function seed_dipersal()
	# wind
	# pollinator dispersed
end

function seed_establishment()
end
end
