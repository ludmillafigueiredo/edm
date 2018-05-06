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
			fmeantemp::Array{Float64,1}
			ftempsd::Array{Float64,1}
			fmeanprec::Array{Float64,1}
			fprecsd::Array{Float64,1}
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
		function landscape_init(landinit::LandPars)

			landscape = WorldCell[]

			for frag in collect(1:landinit.nfrags)

				fragment = WorldCell[]

				fragt = rand(Normal(landinit.fmeantemp[frag],landinit.ftempsd[frag]),1)[1] + tK
				fragp = rand(Normal(landinit.fmeanprec[frag],landinit.fprecsd[frag]),1)[1]
				for y in collect(1:landinit.fylength[frag]), x in 1:(landinit.fxlength[frag])
					newcell = WorldCell(true,
										fragt,
										fragp,
										Dict())
					push!(fragment,newcell)
				end
				fragment = reshape(fragment,landinit.fxlength[frag],landinit.fylength[frag])
				# add the fragment to the landscape structure
				if frag == 1
					landscape = fragment #when empty, landscape cant cat with frag
				else
					landscape = cat(3,landscape, frag)
				end
			end
			 # reshaping is easier then going through every index of a 3D landscap, creating a WorldCell cell there and parameterizing it. x in inner loop matches reshape order
			#TODO check reshaping for landscape with frags of different sizes
			return landscape
		end

		"""
		arealoss!()

		"""
		function arealoss!(landscape,settings::Dict{String,Any})

		end

		"""
		climate!(landscape,t)
		Update temperature and precipitation values according to the weekly input data (weekly means and ).
		"""
		function climate!(landscape::Array{Setworld.WorldCell,3}, t)
			#TODO Unity test
			for frag in collect(1:size(landscape[:,:,3])) #every fragment
				for cell in eachindex(landscape[:,:,frag])
					# weekly update temperature
					landscape[cell,frag].temp = rand(Normal(disturbance[t,1],LandRef[t,2]),1)[1] + tK
					# weekly update precipitation
					landscape[cell,frag].precpt = rand(Normal(disturbance[t,3],LandRef[t,4]),1)[1]
				end
			end

		end

		"""
		fragmentation!()
		"""
		function fragmentation!()
		end

	end
