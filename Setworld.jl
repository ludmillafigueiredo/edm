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

		export LandPars, WorldCell, landscape_init, climate!

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

		#Types
		mutable struct WorldCell
			avail::Bool
			temp::Float64
			precpt::Float64
			neighs::Dict{String,Float64}
		end

		mutable struct PollCell
			floralres::Dict #ind => amount of floral resource projected in the cell
		end

		#WorldCell() = WorldCell(false, 0.0, 0.0, Dict())

		# Functions
		function landscape_init(landpars::LandPars)

			landscape = WorldCell[]

			for frag in collect(1:landpars.nfrags)

				fragment = WorldCell[]

				fragt = rand(Normal(landpars.meantempts[1],landpars.sdtempts[1]),1)[1] + tK
				fragp = rand(Normal(landpars.meanprects[1],landpars.sdprects[1]]),1)[1]
				for y in collect(1:landpars.fylength[frag]), x in 1:(landpars.fxlength[frag])
					newcell = WorldCell(true,
										fragt,
										fragp,
										Dict())
					push!(fragment,newcell)
				end
				fragment = reshape(fragment,landpars.fxlength[frag],landpars.fylength[frag])
				# add the fragment to the landscape structure
				if frag == 1
					landscape = fragment #when empty, landscape cant cat with frag
				else
					landscape = cat(3,landscape, frag)
				end
			end
			 # reshaping is easier then going through every index of a 3D landscape, creating a WorldCell cell there and parameterizing it. x in inner loop matches reshape order
			#TODO check reshaping for landscape with frags of different sizes
			return landscape
		end

		"""
		arealoss!()

		"""
		function arealoss!(landscape,settings::Dict{String,Any})

		end

		"""
		updatetemp!(landscape,t)
		Update temperature and precipitation values according to the weekly input data (weekly means and ).
		"""
		function updatetemp!(landscape::Array{Setworld.WorldCell,3}, t, landpars::LandPars)
			#TODO Unity test
			for frag landscape[:,:,3] #every fragment
				for cell in eachindex(landscape[:,:,frag])
					# weekly update temperature
					landscape[cell,frag].temp = rand(Normal(landpars.meantempts[1],landpars.sdtempts[1]),1)[1] + tK
					# weekly update precipitation
					landscape[cell,frag].precpt = rand(Normal(landpars.meanprects[1],landpars.sdprects[1]])[1]
				end
			end

		end

		"""
		fragmentation!()
		"""
		function fragmentation!(landscape::Array{Setworld.WorldCell,3}, t, settings::Dict{String,Any})
		end

	end
