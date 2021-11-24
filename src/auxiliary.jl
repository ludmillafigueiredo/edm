# Initial trait values is read from an input file and stored for reference in `sppref::SppRef`.
mutable struct SppRef
    sp::Array{String, 1}
    clonality::Dict{String,Bool}
    kernel::Dict{String,String}
    pollen_vector::Dict{String,String}
    self_failoutcross::Dict{String,Bool}
    self_proba::Dict{String,Float64}
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
    updateK!(landscape, settings, t, tdist)
Update the carrying capacity of the landscape (`K`) and of each gridcell (`C_K`).
It is called during initialization (`t` = 1)  and is called again if landscape is disturbed (at `tdist`; `criticalts` keeps track of the updates and outputs the new values).
"""
function updateK!(K::Float64, landscape::BitArray{2}, settings::Dict{String,Any}, t::Int64, land_pars::LandPars)

    criticalts = Array{Int64,N} where N
    # check on which timesteps to write land dims (ifelse() does not work)y
    if land_pars.disturbance == nothing
        criticalts = [1]
    else
        criticalts = append!([1],
		     unique(sort(vcat(broadcast(-, land_pars.disturbance.td, 1), broadcast(+, land_pars.disturbance.td, 1), land_pars.disturbance.td))))
    end

    if t in criticalts
        # output message
        if t == 1
            open(joinpath(settings["outputat"],settings["simID"],"landlog.txt"),"w") do sim
                println(sim, "week\tK\tC_K\thabitatarea\ttotalarea")
            end
        end

        # habitat area
        habitatarea = length(findall(x -> x == true, landscape))  #number of habitat grid cels ggrid cells
        totalarea = prod(size(landscape)) # total number of grid cell. Habitat or not

	K = C_K*habitatarea
        open(joinpath(settings["outputat"],settings["simID"],"landlog.txt"),"a") do sim
            writedlm(sim,[t K C_K habitatarea totalarea])
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
    for sp in SPP_REF.sp

        mean_opt = SPP_REF.temp_opt[sp]
        std_tol = SPP_REF.temp_tol[sp]

        absolute_fitness = max_fitness*exp(-((mean_annual-mean_opt)^2)/(2*(std_tol^2)))
        SPP_REF.fitness[sp] = absolute_fitness
    end
    open(joinpath(settings["outputat"], settings["simID"], "spp_fitness.csv"),"a") do fitnessfile
    	writedlm(fitnessfile, hcat(t, sp, get(SPP_REF.fitness, sp, "NA")) for sp in collect(keys(SPP_REF.fitness)))
    end
end

"""
    setenv!(landscape,t)
Update temperature and precipitation values according to the weekly input data (weekly means and ).
"""
function setenv!(t::Int64, temp_ts)

    T = temp_ts[t, :meantemp][1] + T_K
    
    if rem(t, 52) == 1
	mean_annual = mean(temp_ts[t:(t+51), :meantemp] .+ T_K)
	return T, mean_annual
    else
	return T
    end
    
end

"""
get_dest(loc, dist, theta)
For currrent location `loc`, calculate new location `(xdest, ydest)` based on the distance of dispersal `dist` and angle `theta`.
"""
function get_dest(loc::Tuple, dist::Float64, theta::Float64)

    xdest = loc[1]+round(Int64, dist*cos(theta), RoundNearestTiesAway)
    ydest = loc[2]+round(Int64, dist*sin(theta), RoundNearestTiesAway)
    suitable = checkbounds(Bool, landscape, xdest, ydest) && landscape[xdest, ydest] == true
	    
    newloc = (dest = (xdest, ydest), suitable = suitable)
	    
    return newloc
	    
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
    Bm = B0_MORT*(sum(values(plant.mass))^(-1/4))*exp(-A_E/(BOLTZ*T))
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


function grow_allocate!(plant, b0grow, flowering_ids)

    B_grow = b0grow*((sum(values(plant.mass))-plant.mass["repr"])^(-1/4))*exp(-A_E/(BOLTZ*T)) # only vegetative biomass fuels growth

    new_mass = B_grow*((2*plant.compartsize + plant.compartsize^(3/4))-(sum(values(plant.mass))-plant.mass["repr"]))

    if plant.id in flowering_ids
       plant.mass["repr"] += new_mass
    else
       map(x -> plant.mass[x] += (1/3)*new_mass,
       	   ["leaves", "stem", "root"])
    end

end

# function to vectorize age increase
function age!(plant::Plant)
    plant.age += 1
end

function sort_die!(sp::String, sppcell_fitness::Dict{String,Float64}, plants_cell::Array{Plant,1}, dying_stage::String, plants::Array{Plant,1}, settings::Dict{String,Any}, t::Int64)

    C_K_sp = C_K * (sppcell_fitness[sp]/sum(collect(values(sppcell_fitness))))

    plantscell_sp = filter(x -> x.sp == sp, plants_cell)

    if sum(vcat(map(x->sum(values(x.mass))-x.mass["root"], plantscell_sp),NOT_0)) > 0

        dying_plants = filter(x -> x.stage == dying_stage, plantscell_sp)
	# order inds so smaller can be killed first with pop!
        dying_sorted = sort(dying_plants, by = x -> sum(values(x.mass)), rev = true)

	dying = Plant[]
			   
	while (sum(vcat(map(x->sum(values(x.mass))-x.mass["root"], plantscell_sp),NOT_0)) > C_K_sp
	       && length(dying_sorted) > 0)

	    # kill smallest
	    pop!(dying_sorted) |> x -> (push!(dying, x);
	                                filter!(y -> y.id != x.id, plantscell_sp))
	end

        dying_idxs = findall(x -> x.id in getfield.(dying, :id), plants)
        deleteat!(plants, dying_idxs)
        # check point life history events
        open(joinpath(settings["outputat"],settings["simID"],"events.csv"),"a") do sim
            for i in 1:length(dying)
                writedlm(sim, hcat(t, "death-dep", dying_stage, mean(getfield.(dying, :age)), 1))
            end
        end
    end			
end
