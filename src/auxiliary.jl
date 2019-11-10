# Julia Packages
using ArgParse
using Distributions
using CSV
using DataFrames
using DataValues
using RCall
using Random
using DelimitedFiles
using Dates

# Initial trait values is read from an input file and stored for reference in `sppref::SppRef`.
mutable struct SppRef
    sp_id::Array{String, 1}
    clonality::Dict{String,Bool}
    kernel::Dict{String,String}
    seedmass::Dict{String,Float64}
    compartsize::Dict{String,Float64}
    span_min::Dict{String,Float64}
    span_max::Dict{String,Float64}
    firstflower_min::Dict{String,Float64}
    firstflower_max::Dict{String,Float64}
    floron::Dict{String,Float64}
    floroff::Dict{String,Float64}
    seednumber_min::Dict{String,Float64}
    seednumber_max::Dict{String,Float64}
    seedon::Dict{String,Float64}
    seedoff::Dict{String,Float64}
    bankduration_min::Dict{String,Float64}
    bankduration_max::Dict{String,Float64}
    b0grow::Dict{String,Float64}
    b0germ::Dict{String,Float64}
    b0mort::Dict{String,Float64}
    temp_opt::Dict{String,Float64}
    temp_tol::Dict{String,Float64}
    fitness::Dict{String, Float64}
end

SppRef() = SppRef(String[],
	          Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict())

# Minimal and maximal trait values, which control microevolution, are stored in `traitranges::TraitRanges`
mutable struct TraitRanges
    compartsize::Dict{String,Array{Float64,1}}
    span::Dict{String,Array{Int64,1}}
    firstflower::Dict{String,Array{Int64,1}}
    floron::Dict{String,Array{Int64,1}}
    floroff::Dict{String,Array{Int64,1}}
    seednumber::Dict{String,Array{Int64,1}}
    seedon::Dict{String,Array{Int64,1}}
    seedoff::Dict{String,Array{Int64,1}}
    bankduration::Dict{String,Array{Int64,1}}
end

TraitRanges() = TraitRanges(Dict(),
		            Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict())

function gridsizes(realxlen::Array{Float64,1}, realylen::Array{Float64,1})
  xlength = round(Int64,((realxlen.*100)./100), RoundNearestTiesAway)
  ylength = round(Int64,((realylen.*100)./100), RoundNearestTiesAway)
  return xlength, ylength
end

"""
lengthtocell(d)
Converts distance values from m to cell size.
"""
function lengthtocell(d::Float64)
  celldist = round(Int64,d, RoundNearestTiesAway)
  return celldist
end

"""
areatocell(area)
Calculate the side length of a square grid of a fragment of `area` m².
1 m² = 10000 cm², 3 cm of cell side size
# 2 methods:
- single fragment
- multiple fragments
"""
function areatocell(area::Float64)
  side = round(Int64,((sqrt(area*10000))/100), RoundNearestTiesAway)
  return side
end

function areatocell(area::Array{T,1} where {T<:Number})
  side = round.(Int64,(sqrt.(area*10000))/100, RoundNearestTiesAway)
  return side
end

"""
    updateK!(landavail, settings, t, tdist)
Update the carrying capacity of the landscape (`K`) and of each gridcell (`cK`).
It is called during initialization (`t` = 1)  and is called again if landscape is disturbed (at `tdist`; `criticalts` keeps track of the updates and outputs the new values).
"""
function updateK!(K::Float64, landavail::BitArray{2}, settings::Dict{String,Any}, t::Int64, tdist::Any)

    criticalts = Array{Int64,N} where N
    # check on which timesteps to write land dims (ifelse() does not work)
    if tdist == nothing
        criticalts = [1]
    else
        criticalts = vcat(1, (tdist-1), tdist, (tdist+1))
    end

    if t in criticalts
        # output message
        if t == 1
            open(joinpath(settings["outputat"],settings["simID"],"landlog.txt"),"w") do sim
                println(sim, "week\tK\tcK\thabitatarea\ttotalarea")
            end
        end

        # habitat area
        habitatarea = length(findall(x -> x == true, landavail))  #number of habitat grid cels ggrid cells
        totalarea = prod(size(landavail)) # total number of grid cell. Habitat or not

	K = cK*habitatarea
        open(joinpath(settings["outputat"],settings["simID"],"landlog.txt"),"a") do sim
            writedlm(sim,[t K cK habitatarea totalarea])
        end
    end
    return K
end

"""
    updatefitness!(sppref::SppRef, mean_opt::Float64, std_tol::Float64, mean_annual::Float64, max_fitness::Float64)

Calculate the fitness value according to a Gauss function:

         f(x) = a*exp(-((x-b)²)/(2*c²)),

and store the value in a dictionnary holding the species performance at a given time step.
Fitness is updated every begining of the year, with mean temperature for the year `mean_annual` calculated before hand.

# Arguments
- `mean_annual`: mean annual temperature for the following year. `x` in the Gaussian function takes its value.
- `maximal fitness::Float64=1.0`: parameter `a`, the height of the curve's peak.
- `mean_opt::Float64`: parameter `b` is the position of the center of the peak.
- `std_tol::Float64`: parameter `c` is the standard deviation,
"""
function updatefitness!(mean_annual::Float64, max_fitness::Float64, t::Int64, settings::Dict{String, Any})
    for sp in SPP_REFERENCE.sp_id

        mean_opt = SPP_REFERENCE.temp_opt[sp]
        std_tol = SPP_REFERENCE.temp_tol[sp]

        absolute_fitness = max_fitness*exp(-((mean_annual-mean_opt)^2)/(2*(std_tol^2)))
        SPP_REFERENCE.fitness[sp] = absolute_fitness
    end
    open(joinpath(settings["outputat"], settings["simID"], "spp_fitness.csv"),"a") do fitnessfile
    	writedlm(fitnessfile, hcat(t, sp, get(SPP_REFERENCE.fitness, sp, "NA")) for sp in collect(keys(SPP_REFERENCE.fitness)))
    end
end

"""
    updateenv!(landscape,t)
Update temperature and precipitation values according to the weekly input data (weekly means and ).
"""
function setenv!(t::Int64, landpars::LandPars)

    T = landpars.meantempts[t] + tK
    if rem(t, 52) == 1
	mean_annual = mean(landpars.meantempts[t:(t+51)] + tK)
	#unity test
	println("Temperature for week $t: $T")
	println("Mean for the year of week $t: $mean_annual")

	return T, mean_annual
    else
	return T
    end
end

function setenv!(t::Int64, landpars::NeutralLandPars)

    T = landpars.meantempts[t] + tK
    if rem(t, 52) == 1
	mean_annual = mean(broadcast(+, tK, landpars.meantempts[t:(t+51)]))
	#unity test
	println("Temperature for week $t: $T")
	println("Mean for the year of week $t: $mean_annual")

	return T, mean_annual
    else
	return T
    end
end

"""
get_dest(loc, dist, theta)
For currrent location `loc`, calculate new location `(xdest, ydest)` based on the distance of dispersal `dist` and angle `theta`.
"""
function get_dest(loc::NamedTuple{(:idx, :loc),Tuple{Int64,Tuple{Int64,Int64}}},
	 	  dist::Float64,
		  theta::Float64)

    xdest = loc.loc[1]+round(Int64, dist*cos(theta), RoundNearestTiesAway)
    ydest = loc.loc[2]+round(Int64, dist*sin(theta), RoundNearestTiesAway)
    suitable = checkbounds(Bool, landavail, xdest, ydest) && landavail[xdest, ydest] == true
	    
    newlocs = (idx = loc.idx, dest = (xdest, ydest), suitable = suitable)
	    
    return newlocs
	    
end

"""
"""
function vectorized_seedproc(process::String, processing_sp::Array{Plant,1}, B::Float64)

	prob = 1-exp(-B)

	# unit test
	if prob < 0
	    error("$process probability < 0")
	elseif prob > 1
	    error("$process probability > 1")
	end

	n_procs = rand(Distributions.Binomial(length(processing_sp), prob))[1]

	ids_procs = sample(getfield.(processing_sp, :id), n_procs) 

	return ids_procs # return ids and not indexes because germination uses ids,
	       		 # and mortality uses indexes
end

"""
mort_prob(plant)
Calculate mortality probability for a single plant
"""
function mort_prob(plant::Plant, T)
    Bm = B0_MORT*(sum(values(plant.mass))^(-1/4))*exp(-aE/(Boltz*T))
    mort_prob = 1-exp(-Bm)
    return mort_prob
end

"""
survival(dying)
Calculate the 
"""
function survival(dying::Array{Plant, 1}, T)
    deaths = map(x -> (id = x.id,
                       death = rand(Bernoulli(mort_prob(x, T)))), dying)
    dead_ids = getfield.(filter(x -> x.death == 1, deaths), :id)
    living_ids = getfield.(filter(x -> x.death == 0, deaths), :id)
    return dead_ids, living_ids 
end