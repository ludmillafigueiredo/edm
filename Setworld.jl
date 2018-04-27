"""
This module contains the type of the cell and functions for setting up initial environmental conditions and changing it when necessary.

	# WorldCell type
	This type has fields designed to store all relevant environmental conditions and plants steming points. Other types of organisms are sotred in arrays and have a location, but do not need to "point" to a landscape cell.

		# landscape_init()
		landscape_init(fragmentation, xlength,ylength, n_frags, meantemp, tempsd, meanprec, precsd)

		landscape_init() creates a multidimensional array composed of habitat grid cells (Array{WorldCell}(xlength, ylength, n_frags)). If -fragmentation` is `false`, then one single grid of `xlength` x `ylength` cells is created. If `fragmentation` is `true`, `n_frags`
		# Arguments
		`fragmentation::Bool` statement about fragmentation of the landscape
		`xlength::Int64` Maximal x length of a fragment or total x length of an unfragmented landscape
		`ylength::Int64` Maximal y length of a fragment or total y length of an unfragmented landscape
		`n_frags::Int64` Number of fragments
		`meantemp::Float64` and `tempsd::Float64` mean temperature of a patch and standard deviation
		`meanprec::Float64` and `precsd::Float64` mean precipitation  of a patch and standard-deviation
		"""
		module Setworld

		using Distributions

		export Simpars, Landpars, WorldCell, landscape_init, climate!

		const tK = 273.15 # °C to K converter

		#Simulation parameters storage:
		mutable struct Landpars #TODO put all landscape.in values in here
			fxlength::Array{Int64,N} where N
			fylength::Array{Int64,N} where N
			fmeantemp::Array{Float64,N} where N
			ftempsd::Array{Float64,N} where N
			fmeanprec::Array{Float64,N} where N
			fprecsd::Array{Float64,N} where N
			nfrags::Int64
		end

		#Types
		mutable struct WorldCell
			avail::Bool
			#resourcesFONs::Dict #separate Dict matrix
			#pollinationFONs::Dict  separate Dict matrix
			temp::Float64
			precpt::Float64
			neighs::Dict{String,Float64} #
			#stem::Bool #steming point of plant
			#WorldCell() = new() #if this is included, and WorldCell object must be initialized as WordCell() and then completed
		end

		mutable struct PollCell
			floralres::Dict #ind => amount of floral resource projected in the cell
		end

		#WorldCell() = WorldCell(false, 0.0, 0.0, Dict())

		# Functions
		function landscape_init(simpars)

			landscape = WorldCell[]

			for frag in 1:simpars.nfrags
				for y in 1:simpars.fylength[frag], x in 1:simpars.fxlength[frag]
					newcell = WorldCell(true,
					rand(Normal(simpars.fmeantemp[frag],simpars.ftempsd[frag]),1)[1] + tK,
					rand(Normal(simpars.fmeanprec[frag],simpars.fprecsd[frag]),1)[1],
					Dict())
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
		climate!(landscape,t)
		Update temperature and precipitation values according to the weekly input data (weekly means and ).
		"""
		function climate!(landscape::Array{Setworld.WorldCell,3}, t)
			#TODO Unity test
			for frag in 1:size(landscape[:,:,3]) #every fragment
				for cell in eachindex(landscape[:,:,frag])
					# weekly update temperature
					landscape[cell,frag].temp = rand(Normal(LandRef[t,1],LandRef[t,2]),1)[1] + tK
					# weekly update precipitation
					landscape[cell,frag].precpt = rand(Normal(LandRef[t,3],LandRef[t,4]),1)[1]
				end
			end

		end
	end
