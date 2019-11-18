
"""
    init_K(landavail, settings, t)
Calculate the carrying capacity of the landscape (`K`) and of each gridcell (`C_K`) at initialization.
`K` is used to initialize the species with abundances corresponding to a niche partitioning model.
"""
function init_K(landavail::BitArray{2}, settings::Dict{String,Any}, t::Int64)
    habitatarea = length(findall(x -> x == true, landavail))
    totalarea = prod(size(landavail))
    K = C_K*habitatarea
    return K
end

"""
    init_landscape()
Create the initial landscape structure.
"""
function init_landscape(landpars::LandPars)

    landscape = Array{Dict{String,Float64},2} #if created only inside the loop, remains a local variable

    for p in collect(1:landpars.npatches)

	patch = fill(Dict{String,Float64}(), landpars.plength[p],landpars.plength[p])

	if p == 1
	    landscape = patch #when empty, landscape cant cat with frag
	else
	    landscape = cat(3,landscape, patch)
	end
    end

    landavail = fill(true,size(landscape))

    return landscape, landavail
end

function init_landscape(landpars::NeutralLandPars)

    # convert matrix to BitArray (smaller than Bool)
    landavail = Bool.(landpars.initialland)

    # create landscape with same dimensions
    landscape = fill(Dict{String, Float64}(),
		     size(landavail))

    return landscape, landavail

end

"""
    updatefitness!(SPP_REF, mean_opt::Float64, std_tol::Float64, mean_annual::Float64, max_fitness::Float64)

Calculate the fitness value according to a Gauss function:

         f(x) = a*exp(-((x-b)²)/(2*c²)),

and store the value in a dictionnary holding the species performance at a given time step.
It runs at initialization, to set up niche partitioning. It used the  with mean temperature for the year `mean_annual` calculated before hand

# Arguments
- `mean_annual`: mean annual temperature for the following year. `x` in the Gaussian function takes its value.
- `maximal fitness::Float64=1.0`: parameter `a`, the height of the curve's peak.
- `mean_opt::Float64`: parameter `b` is the position of the center of the peak.
- `std_tol::Float64`: parameter `c` is the standard deviation,
"""
function init_fitness!(SPP_REF::SppRef, mean_annual::Float64, max_fitness::Float64)
    for sp in SPP_REF.sp_id
        mean_opt = SPP_REF.temp_opt[sp]
        std_tol = SPP_REF.temp_tol[sp]

        absolute_fitness = max_fitness*exp(-((mean_annual-mean_opt)^2)/(2*(std_tol^2)))
        SPP_REF.fitness[sp] = absolute_fitness
    end
end


"""
    init_plants(landavail, SPP_REF,id_counter)

Initialize the organisms with trait values stored in `SPP_REF` and distributes them in the suitable grid-cells in the landscape `landavail`.
Store the individuals in the `plants` array, which holds all plants simulated at any given time.
"""
function init_plants(landavail::BitArray{N} where N, SPP_REF::SppRef, id_counter::Int, settings::Dict{String, Any}, K::Float64)

    plants = Plant[]

    for s in SPP_REF.sp_id

        # Niche partitioning: Upon initialization, each species total biomass equals K*fitness_relative, where fitness_relative is the species fitness values relative to the sum of others
        # The initial abundance is number of medium-sized individuals (50% of maximal biomass) that would sum up to the biomass.
        sp_abund = Int(round(((SPP_REF.fitness[s]/sum(collect(values(SPP_REF.fitness))))*K)/(0.5*(2*SPP_REF.compartsize[s]+SPP_REF.compartsize[s])), RoundUp))

	open(joinpath(settings["outputat"], string(settings["simID"], "initialabundances.txt")),"a") do sim
	    println(sim, "Initial abundance of $s: $sp_abund")
        end

	XYs = hcat(rand(1:size(landavail,1), sp_abund),
                    rand(1:size(landavail,2), sp_abund))

	for i in 1:sp_abund

	    id_counter += 1 # update individual counter
	    minvalue = 1e-7 # Distribution.Uniform requires max > min

	    newplant = Plant(string(id_counter, base=16),
			      rand(["a" "j" "s"]),
			      (XYs[i,1],XYs[i,2]),
			      s,
			      SPP_REF.kernel[s],
			      SPP_REF.clonality[s],
			      SPP_REF.compartsize[s], #compartsize
			      Int(round(rand(Distributions.Uniform(SPP_REF.span_min[s], SPP_REF.span_max[s] + minvalue),1)[1], RoundUp)),
			      Int(round(rand(Distributions.Uniform(SPP_REF.firstflower_min[s], SPP_REF.firstflower_max[s] + minvalue),1)[1], RoundUp)),
			      SPP_REF.floron[s],
			      SPP_REF.floroff[s],
			      Int(round(rand(Distributions.Uniform(SPP_REF.seednumber_min[s],SPP_REF.seednumber_max[s] + minvalue),1)[1], RoundUp)),
			      SPP_REF.seedon[s],
			      SPP_REF.seedoff[s],
			      Int(round(rand(Distributions.Uniform(SPP_REF.bankduration_min[s],SPP_REF.bankduration_max[s] + minvalue),1)[1], RoundUp)),
			      0, #age
			      Dict("leaves" => 0.0, "stem" => 0.0, "repr" => 0.0, "root" => 0.0),
			      false)

            ## initial biomass
	    if newplant.stage == "s"
		newplant.mass["root"] = SPP_REF.seedmass[newplant.sp]
		newplant.age = 1
	    elseif newplant.stage == "j"
		newplant.mass["root"] = SPP_REF.seedmass[newplant.sp]
		newplant.age = newplant.seedon + 4
	    elseif newplant.stage in ["a"]
		newplant.mass["leaves"] = newplant.compartsize^(3/4) * 0.75
		newplant.mass["stem"] = newplant.compartsize * 0.75
		newplant.mass["root"] = newplant.compartsize * 0.75
		newplant.age = newplant.firstflower
	    else
		error("Check individual stages.")
	    end

	    push!(plants, newplant)

	end

    open(joinpath(settings["outputat"], string(settings["simID"], "initialabundances.txt")),"a") do sim
        println(sim, "Seeds of $s: $(length(findall(x -> x.sp == s && x.stage == "s", plants)))")
	println(sim, "Juveniles of $s: $(length(findall(x -> x.sp == s && x.stage == "j", plants)))")
	println(sim, "Adults of $s: $(length(findall(x -> x.sp == s && x.stage == "a", plants)))")
    end
	
    end
    return plants, id_counter
end

settings = parse_commandline()
println(keys(settings))
Random.seed!(settings["rseed"])

id_counter = 0
management_counter = 0

tdist = set_tdist(settings)
landpars = read_landpars(settings)

const SPP_REF = read_sppinput(settings)
const TRAIT_RANGES = define_traitranges(settings)
interaction, scen, remaining = read_insects(settings)
global mylandscape, landavail = init_landscape(landpars)
K = init_K(landavail, settings, 1)
T, mean_annual = setenv!(1, landpars)
init_fitness!(SPP_REF, mean_annual, 1.0)
plants, id_counter = init_plants(landavail, SPP_REF, id_counter, settings, K)

#check-points
println("Landscape initialized: type $(typeof(mylandscape))")
println("Landscape is object of type $(typeof(landpars))")
println("Plants initialized: type $(typeof(plants))")
println("Starting simulation")