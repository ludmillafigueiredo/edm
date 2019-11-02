
"""
    init_K(landavail, settings, t)
Calculate the carrying capacity of the landscape (`K`) and of each gridcell (`cK`) at initialization.
`K` is used to initialize the species with abundances corresponding to a niche partitioning model.
"""
function init_K(landavail::BitArray{2}, settings::Dict{String,Any}, t::Int64)
    habitatarea = length(findall(x -> x == true, landavail))
    totalarea = prod(size(landavail))
    K = cK*habitatarea
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
    updatefitness!(sppref, mean_opt::Float64, std_tol::Float64, mean_annual::Float64, max_fitness::Float64)

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
function init_fitness!(sppref::SppRef, mean_annual::Float64, max_fitness::Float64)
    for sp in sppref.sp_id
        mean_opt = sppref.temp_opt[sp]
        std_tol = sppref.temp_tol[sp]

        absolute_fitness = max_fitness*exp(-((mean_annual-mean_opt)^2)/(2*(std_tol^2)))
        sppref.fitness[sp] = absolute_fitness
    end
end


"""
    init_plants(landavail, sppref,id_counter)

Initialize the organisms with trait values stored in `sppref` and distributes them in the suitable grid-cells in the landscape `landavail`.
Store the individuals in the `plants` array, which holds all plants simulated at any given time.
"""
function init_plants(landavail::BitArray{N} where N, sppref::SppRef, id_counter::Int, settings::Dict{String, Any}, K::Float64)

    plants = Plant[]

    for s in sppref.sp_id

        # Niche partitioning: Upon initialization, each species total biomass equals K*fitness_relative, where fitness_relative is the species fitness values relative to the sum of others
        # The initial abundance is number of medium-sized individuals (50% of maximal biomass) that would sum up to the biomass.
        sp_abund = Int(round(((sppref.fitness[s]/sum(collect(values(sppref.fitness))))*K)/(0.5*(2*sppref.compartsize[s]+sppref.compartsize[s])), RoundUp))

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
			      sppref.kernel[s],
			      sppref.clonality[s],
			      sppref.compartsize[s], #compartsize
			      Int(round(rand(Distributions.Uniform(sppref.span_min[s], sppref.span_max[s] + minvalue),1)[1], RoundUp)),
			      Int(round(rand(Distributions.Uniform(sppref.firstflower_min[s], sppref.firstflower_max[s] + minvalue),1)[1], RoundUp)),
			      sppref.floron[s],
			      sppref.floroff[s],
			      Int(round(rand(Distributions.Uniform(sppref.seednumber_min[s],sppref.seednumber_max[s] + minvalue),1)[1], RoundUp)),
			      sppref.seedon[s],
			      sppref.seedoff[s],
			      Int(round(rand(Distributions.Uniform(sppref.bankduration_min[s],sppref.bankduration_max[s] + minvalue),1)[1], RoundUp)),
			      0, #age
			      Dict("leaves" => 0.0, "stem" => 0.0, "repr" => 0.0, "root" => 0.0),
			      false)

            ## initial biomass
	    if newplant.stage == "s"
		newplant.mass["root"] = sppref.seedmass[newplant.sp]
		newplant.age = newplant.seedon + 1
	    elseif newplant.stage == "j"
		newplant.mass["root"] = sppref.seedmass[newplant.sp]
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
    end
    return plants, id_counter
end

settings = parse_commandline()
println(keys(settings))
Random.seed!(settings["rseed"])

tdist = nothing
if settings["disturbtype"] != "none"
    tdist = loaddisturbance(settings)
end
repr(tdist)

id_counter = 0
management_counter = 0

landpars = read_landpars(settings)
sppref = read_spinput(settings)
traitranges = define_traitranges(settings)
interaction, scen, remaining = read_insects(settings)
global mylandscape, landavail = init_landscape(landpars)
K = init_K(landavail, settings, 1)
T, mean_annual = setenv!(1, landpars)
init_fitness!(sppref, mean_annual, 1.0)
plants, id_counter = init_plants(landavail, sppref, id_counter, settings, K)

#check-points
println("Landscape initialized: type $(typeof(mylandscape))")
println("Landscape is object of type $(typeof(landpars))")
println("Plants initialized: type $(typeof(plants))")
println("Starting simulation")