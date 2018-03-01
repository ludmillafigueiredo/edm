cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/")

module plants

# Load packages
using DataFrames #depreciated readtable fct
using CSV

# Public types, functions and export

#include

mutable struct Plant
	id::String
	sp::String
	location::Tuple # fixed, cant change
	PFT::String
	genotype::Array{String}
	 #TODO own type? traits depend on them
	#TODO traits: traits depend on genotype
	vegmass::Float64
	floral_mass::Float64
	pollen_disp_synd::String #wind dispersal or pollinator
	disp_par::Tuple #kernel parameters for wind dipersal. Pollinator dispersal is done by popllinator dispersal and survival
	poll_offer::String #temp: should match the pollinator's preference in "resouce_preference". Should it be a float, for trait matching?
	#pollmatchtrait #TODO for trait matchingS.
	#TODO evolutionary change should increase prob of trait matching (Poisot Beyond)
	c_FON::Float64
	Plant() = new()
end

function read_plantslist(plantfile::String)
	plant_info = CSV.read(plantfile, header = true)
	return plant_info
end


global plant_id_counter = Int64(0) #TODO after being initialized, it will have to keep track of all other new individuals

function createPlant(x::Int64,y::Int64)
	# 1. Create sp acording to species data: #TODO use different methods to generate from seeds
	plant_list = read_plantslist(plantfile)

	+(plant_id_counter,1)

	newplant = Plant()
	newplant.id = string(plant, plant_id_counter)
	newplant.sp = rand(GLOBAL_RNG, plant_list[:sps],1)[1] #index because otherwise it is an array
	newplant.PFT = plant_list[plant_list[:sps] .== newplant.sp, :ptf][1]
	#newplant.genotype = ["A";"A"] #TODO own type? traits depend on them
	#TODO traits: traits depend on genotype
	newplant.vegmass = 0.1, #TODO read from file
	newplant.floral_mass = 0, #TODO read from file
	#poll_offer = 'flower',#TODO read from file
	#pollen_disp_synd = 'poll' ,#TODO read from file: wind dispersal or pollinator
	#disp_par = [1 1], #TODO read from file, kernel parameters
	# fixed, cant change
	#newplant.c_FON = 1 #TODO make it dependent on biomass

	#store id somewhere
	# Project the neighborhood
	FONproj(x,y,landscape, Plant) # TODO different method for pollination FON

	return newplant, indplants, plant_id_counter
end

function populate(xlength,ylength, agent::Plant)
	#TODO create method for pollinators and herbirvores? They move! Not a fixed(tuple) location field? 
	# spread individuals in the FRAGMENT (not landscape)

	return location
end

#TODO store all plants somewhere or mk it easier to retrieve
