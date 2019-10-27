"""
This module contains the data structures and functions to simulate plants, pollination regimes, and life cycle processes:
- allocate!(): biomass growth and resource allocation
- develop!(): maturation of juveniles that have reached their respective age of first flowering
- mate!(): pollination
- mkoffsrping!(): seed/clone production
- getreleases(): find seeds that can be dispersed
- disperse!(): set new locations for seeds being dispersed
- establish!(): see if seeds manage to germinate in current location
- shedd!(): decrease of reproductive biomass at the end of reproductive season, or winter die-back
- manage!(): decrease of vegetative and reproductive biomass due to mowing

It also contains the data structures and functions for setting up initial environmental conditions, the simulation grid,  and changing it when necessary.
- LandPars and NeutralLandPars: store landscape dimensions, temperature time-series and time(s) of disturbance
- landscape_init(): create the initial simulation grid
- updateenv!(): update weekly temperature and annual mean
- destroyarea!(): make contiguous area of landscape unsuitable for occupation
- fragment!(): make non-contiguous area of landscape unsuitable for occupation
"""
module submodels

using Distributions
using DataFrames
using DataValues
using StatsBase
using DelimitedFiles
using auxfunctions
using entities

include("constants_globalpars.jl")

export initplants, develop!, allocate!, mate!, mkoffspring!, microevolution!, disperse!, germinate, establish!, survive!, shedflower!, winter_dieback!, manage!, destroyorgs!, getreleases, landscape_init, updateenv!, destroyarea!, fragment!

"""
    initplants(landavail, sppref,id_counter)

Initialize the organisms with trait values stored in `sppref` and distributes them in the suitable grid-cells in the landscape `landavail`.
Store the individuals in the `plants` array, which holds all plants simulated at any given time.
"""
function initplants(landavail::BitArray{N} where N, sppref::SppRef, id_counter::Int, settings::Dict{String, Any}, K::Float64)

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
			      sppref.seedmass[s],
			      sppref.compartsize[s], #compartsize
			      Int(round(rand(Distributions.Uniform(sppref.span_min[s], sppref.span_max[s] + minvalue),1)[1], RoundUp)),
			      Int(round(rand(Distributions.Uniform(sppref.firstflower_min[s], sppref.firstflower_max[s] + minvalue),1)[1], RoundUp)),
			      sppref.floron[s],
			      sppref.floroff[s],
			      Int(round(rand(Distributions.Uniform(sppref.seednumber_min[s],sppref.seednumber_max[s] + minvalue),1)[1], RoundUp)),
			      sppref.seedon[s],
			      sppref.seedoff[s],
			      Int(round(rand(Distributions.Uniform(sppref.bankduration_min[s],sppref.bankduration_max[s] + minvalue),1)[1], RoundUp)),
			      sppref.b0grow[s],
			      sppref.b0germ[s],
			      sppref.b0mort[s],
			      0, #age
			      Dict("leaves" => 0.0, "stem" => 0.0, "repr" => 0.0, "root" => 0.0),
			      false)

            ## initial biomass
	    if newplant.stage == "s"
		newplant.mass["root"] = newplant.seedmass
		newplant.age = newplant.seedon + 1
	    elseif newplant.stage == "j"
		newplant.mass["root"] = newplant.seedmass
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

"""
    allocate!(orgs, t, aE, Boltz, setting, sppref, T)
"""
function allocate!(plants::Array{Plant,1}, t::Int64, aE::Float64, Boltz::Float64, settings::Dict{String, Any},sppref::SppRef, T::Float64, biomass_production::Float64, K::Float64, growing_stage::String)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Growth of $growing_stage")
    end

    growing = findall(x->x.stage == growing_stage,plants)

    for o in growing

        b0grow = plants[o].b0grow

	current_vegmass = plants[o].mass["leaves"] + plants[o].mass["stem"] + plants[o].mass["root"]

	B_grow = b0grow*(current_vegmass^(-1/4))*exp(-aE/(Boltz*T)) # only vegetative biomass fuels growth

	new_mass = B_grow*((2*plants[o].compartsize + plants[o].compartsize^(3/4))-current_vegmass)

	# adults of a minimal sie in their reproductive season allocate to reproductive structures, instead of leaves and stem
        # otherwise, growth is equally divided between all vegetative structures
        if (plants[o].stage == "a" &&
            (plants[o].floron <= rem(t,52) < plants[o].floroff) &&
            current_vegmass >= 0.5*(2*plants[o].compartsize+plants[o].compartsize^(3/4)))

	    if haskey(plants[o].mass,"repr")
         	plants[o].mass["repr"] += new_mass
            else
		plants[o].mass["repr"] = new_mass
            end
	else
            plants[o].mass["leaves"] += (1/3)*new_mass
	    plants[o].mass["stem"] += (1/3)*new_mass
	    plants[o].mass["root"] += (1/3)*new_mass
        end
    end
    # unit test
    masserror = findall(x -> sum(collect(values(x.mass))) <= 0, plants)
    if length(masserror) > 0
        println("Plant: plants[masserror[1]]")
	error("Zero or negative values of biomass detected.")
    end
end

"""
    develop!(plants, settings, t)
Juveniles in `plants` older than their age of first flowering at time `t` become adults.
"""
function develop!(plants::Array{submodels.Plant,1}, settings::Dict{String, Any}, t::Int)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Maturation...")
    end

    # find indexes of individuals that are ready to become adults
    juvs = findall(x->x.stage == "j" && x.age >= x.firstflower,plants)

    for j in juvs
	plants[j].stage = "a"
        # check-point
	open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
	    writedlm(sim, hcat(t, "maturation", plants[j].stage, plants[j].age))
	end
    end
end

"""
    mate!(plants, t, setting, scen, tdist, remaining)
Calculate proportion of `plants` that reproduced at time `t`, acording to pollination scenario `scen`, and mark that proportion of the population with the `mated` label.
"""
function mate!(plants::Array{submodels.Plant,1}, t::Int, settings::Dict{String, Any}, scen::String, tdist::Any, remaining)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Pollination ...")
    end

    ready = findall(x-> x.stage == "a" && x.mass["repr"] > x.seedmass, plants)
    pollinated = []
    npoll = 0

    if length(ready) > 0 # check if there is anyone flowering

        # Scenarios were pollination is not species-specific
        # --------------------------------------------------
        if scen in ["indep" "equal"]

	    npoll_default = rand(Distributions.Binomial(Int(ceil(length(ready) * VISITED_DEFAULT)), 0.6))[1]

	    # Determine number of individuals that get pollinated (species is not relevant)
	    if scen == "indep"
		npoll = npoll_default # Fishman & Hadany's proportion of visited flowers
		println("Scenario of INDEP pollination loss")
		# check-point
		open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        	    println(sim, "Number of pollinated: $npoll")
    		end

	    elseif scen == "equal" #all species lose pollination randomly (not species-specific)
		println("Scenario of EQUAL pollination loss")
		if t in tdist                                                          # calculate the amount of loss for the specified times
		    npoll = Int(ceil(npoll_default * remaining[findall(tdist == [t])[1]]))
		else
		    npoll = npoll_default
		end
	    end

	    # unity test
        if npoll < 0
		error("Negative number of plants being pollinated.")
	    end

	    if npoll > 0       # check if any should actually be pollinated
		# who are they
		pollinated = sample(ready, npoll, replace = false, ordered = true)
		# pollinate
		for p in pollinated
		    plants[p].mated = true
		end
	    end
	elseif scen in ["rdm" "spec"] #not yet tested

	    # determine which species will loose pollination
	    if scen == "rdm"
		# randomly pick plant species that will loose pollination (n = nspppoll) at a given timestep and find their number
		# pseudo: rdmly pick species:
		# pseudo: spppoll = unique(getfields(plants, :sp)) |> sample(, nspppoll)
	    elseif scen == "spec" #not yet tested
		# pseudo: from a list of loss pollinators, find the plant species (in the interaction matrix) that will loose pollination at a given timestep
	    end

	    # Species-specific pollination
            # ----------------------------
            # check if any should be pollinated
	    if npoll == 0

	    elseif npoll > 1
		pollinated = sample(ready, npoll, replace = false, ordered = true)
		for p in pollinated
		    plants[p].mated = true
		end
	    else
		error("Negative number of plants being pollinated.")
	    end
	else
	    error("Please chose a pollination scenario \"scen\" in insect.csv:
                   - \"indep\": sexual reproduction happens independently of pollination
                   - \"rmd\": random loss of pollinator species (complementary file should be provided, see model dodumentation)
                   - \"spec\": specific loss of pollinator species (complementary files should be provided, see model documentation)")
	end
    end
end

"""
    mkoffspring!()
After mating happened (marked in `reped`), calculate the amount of offspring each individual produces, both sexually and assexually.
"""
function mkoffspring!(plants::Array{submodels.Plant,1}, t::Int64, settings::Dict{String, Any},sppref::SppRef, id_counter::Int, landavail::BitArray{2}, T::Float64, traitranges::submodels.TraitRanges)
	# check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	println(sim, "Total number of individuals before REPRODUCTION: $length(plants)")
    end

    offspring = Plant[]
    non0sd = 1e-7
    seed_counter = 0

    # Separate sexually and asexually reproducing (mated status can change during the simulation and this would generate more clones)
    ferts = filter(x -> x.mated == true, plants)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
   	 println(sim, "Number of ferts $(length(ferts)) in $(length(plants)) plants")
    end

    asexuals = filter(x -> x.mated == false && x.clonality == true && x.mass["repr"] > x.seedmass, plants)
    
    # Sexuallly produced offspring
    # ----------------------------
    for sp in unique(getfield.(ferts, :sp))

	sowing = findall(x -> x.sp == sp && x.id in getfield.(ferts, :id), plants)
	# check-point
    	open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	    println(sim, "Number of sowing: $(length(sowing))")
        end
	    
	spoffspringcounter = 0 #offspring is not output in the same file as adults and juveniles

	for s in sowing

	    seedmass = plants[s].seedmass
	    offs = div(ALLOC_SEED*plants[s].mass["repr"], seedmass)
	    
	    if offs <= 0
		continue
	    else
		# limit offspring production to the maximal number of seeds the species can produce
		offs > plants[s].seednumber ? offs = plants[s].seednumber : offs

                spoffspringcounter += offs

		# update available biomass for reproduction
		plants[s].mass["repr"] -= (offs * seedmass)

		# unity test
		if  plants[s].mass["repr"] <= 0
		    error("Negative reproductive biomass")
		end

		# get another parent
		conspp = rand(ferts)

                for n in 1:offs

		    id_counter += 1

		    seed = deepcopy(plants[s])
		    seed_counter += 1

		    # unity test
		    for f in fieldnames(typeof(seed))
			if typeof(getfield(seed, f)) in [Int64 Float64]
			    if getfield(seed, f) < 0 || getfield(seed, f) == Inf
				error(f, " has value: ", getfield(seed, f), "Ind: ", seed.id, "sp: ", seed.sp)
			    end
			end
		    end

                    # State variables
		    seed.id = string(id_counter, base = 16)
		    seed.mass = Dict("leaves" => 0.0,
		                     "stem" => 0.0,
				     "root" => seed.seedmass,
				     "repr" => 0.0)
                    seed.stage = "s"
                    seed.age = 0
                    seed.mated = false

                    # Trait microevolution
                    # --------------------
		    #seed.seedmass += rand(Distributions.Normal(0, abs(plants[s].seedmass-conspp.seedmass+non0sd)/6))[1]
			#TODO use fieldnames to loop through  species properties (necessary in other places in the code)
            seed.compartsize = (seed.compartsize+conspp.compartsize)/2+ rand(Distributions.Normal(0, abs(plants[s].compartsize-conspp.compartsize+non0sd)/6))[1]
		    seed.span = Int(round((seed.span+conspp.span)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].span-conspp.span+non0sd)/6))[1], RoundUp))
		    seed.firstflower = Int(round((seed.firstflower+conspp.firstflower)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].firstflower-conspp.firstflower+non0sd)/6))[1], RoundUp))
		    seed.floron = Int(round((seed.floron+conspp.floron)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].floron-conspp.floron+non0sd)/6))[1],RoundUp))
		    seed.floroff = Int(round((seed.floroff+conspp.floroff)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].floroff-conspp.floroff+non0sd)/6))[1],RoundUp))
		    seed.seednumber = Int(round((seed.seednumber+conspp.seednumber)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].seednumber-conspp.seednumber+non0sd)/6))[1], RoundUp))
		    seed.seedon = Int(round((seed.seedon+conspp.seedon)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].seedon-conspp.seedon+non0sd)/6))[1],RoundUp))
		    seed.seedoff = Int(round((seed.seedoff+conspp.seedoff)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].seedoff-conspp.seedoff+non0sd)/6))[1],RoundUp))
		    seed.bankduration = Int(round((seed.bankduration+conspp.bankduration)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].bankduration-conspp.bankduration+non0sd)/6))[1],RoundUp))

		    # constrain microevolution: avoid trait changes that generate negative values (and also values that irrealistically high)

                    #if (seed.seedmass < traitranges.seedmass[seed.sp][1] || seed.seedmass > traitranges.seedmass[seed.sp][end])
		    #   seed.seedmass < traitranges.seedmass[seed.sp][1] ? seed.seedmass = traitranges.seedmass[seed.sp][1] :
		    #   seed.seedmass = traitranges.seedmass[seed.sp][end]
		    #end

		    if (seed.compartsize < traitranges.compartsize[seed.sp][1] || seed.compartsize > traitranges.compartsize[seed.sp][end])
			seed.compartsize < traitranges.compartsize[seed.sp][1] ? seed.compartsize = traitranges.compartsize[seed.sp][1] :
			    seed.compartsize = traitranges.compartsize[seed.sp][end]
		    end

		    if (seed.span < traitranges.span[seed.sp][1] || seed.span > traitranges.span[seed.sp][end])
			seed.span < traitranges.span[seed.sp][1] ? seed.span = traitranges.span[seed.sp][1] :
			    seed.span = traitranges.span[seed.sp][end]
		    end

		    if (seed.firstflower < traitranges.firstflower[seed.sp][1] || seed.firstflower > traitranges.firstflower[seed.sp][end])
			seed.firstflower < traitranges.firstflower[seed.sp][1] ? seed.firstflower = traitranges.firstflower[seed.sp][1] :
			    seed.firstflower = traitranges.firstflower[seed.sp][end]
		    end

		    if (seed.floron < traitranges.floron[seed.sp][1] || seed.floron > traitranges.floron[seed.sp][end])
			seed.floron < traitranges.floron[seed.sp][1] ? seed.floron = traitranges.floron[seed.sp][1] :
			    seed.floron = traitranges.floron[seed.sp][end]
		    end

		    if (seed.floroff < traitranges.floroff[seed.sp][1] || seed.floroff > traitranges.floroff[seed.sp][end])
			seed.floroff < traitranges.floroff[seed.sp][1] ? seed.floroff = traitranges.floroff[seed.sp][1] :
			    seed.floroff = traitranges.floroff[seed.sp][end]
		    end

		    if (seed.seednumber < traitranges.seednumber[seed.sp][1] || seed.seednumber > traitranges.seednumber[seed.sp][end])
			seed.seednumber < traitranges.seednumber[seed.sp][1] ? seed.seednumber = traitranges.seednumber[seed.sp][1] :
			    seed.seednumber = traitranges.seednumber[seed.sp][end]
		    end

		    if (seed.seedon < traitranges.seedon[seed.sp][1] || seed.seedon > traitranges.seedon[seed.sp][end])
			seed.seedon < traitranges.seedon[seed.sp][1] ? seed.seedon = traitranges.seedon[seed.sp][1] :
			    seed.seedon = traitranges.seedon[seed.sp][end]
		    end

		    if (seed.seedoff < traitranges.seedoff[seed.sp][1] || seed.seedoff > traitranges.seedoff[seed.sp][end])
			seed.seedoff < traitranges.seedoff[seed.sp][1] ? seed.seedoff = traitranges.seedoff[seed.sp][1] :
			    seed.seedoff = traitranges.seedoff[seed.sp][end]
		    end

		    if (seed.bankduration < traitranges.bankduration[seed.sp][1] || seed.bankduration > traitranges.bankduration[seed.sp][end])
			seed.bankduration < traitranges.bankduration[seed.sp][1] ? seed.bankduration = traitranges.bankduration[seed.sp][1] :
			    seed.bankduration = traitranges.bankduration[seed.sp][end]
		    end
push!(plants, seed)
end
plants[s].mated = false # after producing seeds in a week, the plant will only do it again in the next week if it gets pollinated again
end
end

open(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv"),"a") do seedfile
    writedlm(seedfile, hcat(t, sp, "s", "sex", spoffspringcounter))
end
end
open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
    writedlm(sim, hcat("Total number of individuals after SEX:", length(plants)))
end

# Asexually produced offspring
# ----------------------------
for sp in unique(getfield.(asexuals, :sp))

    cloning = findall(x -> x.sp == sp && x.id in unique(getfield.(asexuals, :id)) , plants) # find mothers cloning

    # start production counting
    spclonescounter = 0

    for c in cloning

    if rand(Distributions.Binomial(1, 0.5)) == 1
	# get a copy of the mother, which the clones will look like
	clone = deepcopy(plants[c])
	clone.stage = "j" #clones have already germinated
	clone.mass["leaves"] = (plants[c].compartsize*0.1)^(3/4)
	clone.mass["stem"] = plants[c].compartsize*0.1
	clone.mass["root"] = plants[c].compartsize*0.1
	clone.mass["repr"] = 0.0

	id_counter += 1
	clone.id = hex(id_counter)
	push!(plants, clone)
	spclonescounter += 1
    end
    end
    # output clones
    open(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv"),"a") do seedfile
	writedlm(seedfile, hcat(t, sp, "j", "asex", spclonescounter))
    end
end
open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
    writedlm(sim, hcat("Total number of individuals after ASEX:", length(plants)))
end
return id_counter
end


"""
    getreleases(plants, t)
Retrieve indexes of seeds that are still at the mother plant who have entered their sowing season at time `t`.
"""
function getreleases(plants::Array{submodels.Plant,1}, t::Int)
    seedsi = findall(x -> x.stage == "s" && x.age == 0 && x.seedon <= rem(t,52) < x.seedoff, plants)
end

"""
    disperse!(landscape, plants, t, seetings, sppref,)
Seeds are dispersed.

"""
function disperse!(landavail::BitArray{2}, seedsi, plants::Array{submodels.Plant, 1}, t::Int, settings::Dict{String, Any}, sppref::SppRef, landpars::Any, tdist::Any)#submodels.LandPars)}

    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        writedlm(sim, hcat("Dispersing ..."))
    end

    lost = Int64[]
    justdispersed = String[]

    for d in seedsi # only seeds that have been released can disperse

        if plants[d].kernel == "short"
	    dist = auxfunctions.lengthtocell(4*(rand(Distributions.InverseGaussian(µ_short, λ_short),1)[1]))
	elseif plants[d].kernel == "medium"
	    dist = auxfunctions.lengthtocell(1000*(rand(Distributions.InverseGaussian(µ_medium, λ_medium),1)[1]))
	elseif plants[d].kernel == "long"
	    dist = auxfunctions.lengthtocell(rand(Distributions.InverseGaussian(µ_long, λ_long),1)[1])
	elseif plants[d].kernel in ["medium-short", "short-medium"]
	    dist = rand([auxfunctions.lengthtocell(4*(rand(Distributions.InverseGaussian(µ_short, λ_short),1)[1])),
			 auxfunctions.lengthtocell(1000*(rand(Distributions.InverseGaussian(µ_medium, λ_medium),1)[1]))])
	elseif plants[d].kernel in ["medium-long", "long-medium"]
	    dist = rand([auxfunctions.lengthtocell(rand(Distributions.InverseGaussian(µ_long, λ_long),1)[1]),
			 auxfunctions.lengthtocell(1000*(rand(Distributions.InverseGaussian(µ_medium, λ_medium),1)[1]))])
	elseif plants[d].kernel in ["long-short", "short-long"]
	    dist = rand([auxfunctions.lengthtocell(4*(rand(Distributions.InverseGaussian(µ_short, λ_short),1)[1])),
			 auxfunctions.lengthtocell(rand(Distributions.InverseGaussian(µ_long, λ_long),1)[1])])
	elseif plants[d].kernel == "any"
	    dist = rand([auxfunctions.lengthtocell(4*(rand(Distributions.InverseGaussian(µ_short, λ_short),1)[1])),
			auxfunctions.lengthtocell(1000*(rand(Distributions.InverseGaussian(µ_medium, λ_medium),1)[1])),
			auxfunctions.lengthtocell(rand(Distributions.InverseGaussian(µ_long, λ_long),1)[1])])
	else
            # unity test
	    error("Check dispersal kernel input for species $(plants[d].sp).")
	end

	# Find the cell to which it is dispersing
	θ = rand(Distributions.Uniform(0,2),1)[1]*pi
	xdest = plants[d].location[1] + dist*round(Int64, cos(θ), RoundNearestTiesAway)
	ydest = plants[d].location[2] + dist*round(Int64, sin(θ), RoundNearestTiesAway)

	if checkbounds(Bool, landavail, xdest, ydest) && landavail[xdest, ydest] == true  # Check if individual fall inside the habitat area, otherwise, discard it already
                                                                                          # checking the suitability first would make more sense but cant be done if cell is out of bounds
	    plants[d].location = (xdest,ydest)
	    #indexes change when inds get deleted. ID is the only trustworthy way to find them
	    push!(justdispersed, plants[d].id)

	else # if the new location is in an unavailable habitat or outside the landscape, the seed dies

	    push!(lost,d)
	    # check-point of life-history processes
	    open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
		writedlm(sim, hcat(t, "lost in dispersal", plants[d].stage, plants[d].age))
	    end
	end
    end

    # update who actually dispersed
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	println(sim, "Number of dispersing: $(length(justdispersed))")
    end
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	writedlm(sim, hcat("Lost in dispersal:", length(lost)))
    end
    deleteat!(plants,lost)

    return justdispersed
end

"""
"""
function factorized_seedproc(process::String, processing_plants::Array{submodels.Plant,1}, B::Float64, sp::String, plants::Array{submodels.Plant,1})

	prob = 1-exp(-B)

	# unit test
	if prob < 0
	    error("$process probability < 0")
	elseif prob > 1
	    error("$process probability > 1")
	end

	processing_sp = filter(x -> x.sp == sp, processing_plants)
	n_procs = rand(Distributions.Binomial(length(processing_sp), prob))[1]

	ids_procs = sample(getfield.(processing_sp, :id), n_procs)
	idxs_procs = findall(x -> x.id in ids_procs, plants)

	return idxs_procs
end

"""
    establish!
Seed that have already been released (in the current time step, or previously - this is why `seedsi` does not limit who get to establish) and did not die during dispersal can establish in the grid-cell they are in. Germinated seeds mature to juveniles immediately. Seeds that don't germinate stay in the seedbank.

"""
function establish!(justdispersed, plants::Array{submodels.Plant,1}, t::Int, settings::Dict{String, Any}, sppref::SppRef, T::Float64, biomass_production::Float64, K::Float64)

    lost = Int64[]
    Bg = 0
    prob_g = 0

    establishing = filter(x -> (x.id in justdispersed || (x.stage == "s" && x.age >= 1)), plants)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	writedlm(sim, hcat("Seeds trying to ESTABLISH:", length(establishing)))
    end
    
    for sp in unique(getfield.(establishing, :sp))

    	Bg = sppref.b0germ[sp]*(sppref.seedmass[sp]^(-1/4))*exp(-aE/(Boltz*T))

	germs = factorized_seedproc("germination", establishing, Bg, sp, plants)
	
	for g in germs
        
	    plants[g].stage = "j"
	    plants[g].mass["root"] = plants[g].seedmass
	    open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
	        writedlm(sim, hcat(t, "germination", plants[g].stage, plants[g].age))
	    end
        end
    end
end

"""
    survive!(plants, t, cK, K, settings, sppref, landavail, biomass_production, dying_stage)
Plants density-independent mortality is calculated according to the metabolic theory.
Density-dependent mortality is calculated for cells with biomass over the grid-cell carrying capacity `cK`. It kills smaller individuals until the biomass of the grid-cell is above `cK` again.
Density-dependent mortality is not calculated for seeds. Their density-dependent mortality is calculated alongside adults just for convenience.
Both types of mortalities of adults and juveniles are calculated separately (as set by `dying_stage`).
    survive!(plants, t, settings, sppref)
Calculate seed mortality, factorized for all seeds of a same species.
"""
function survive!(plants::Array{submodels.Plant,1}, t::Int, settings::Dict{String, Any}, sppref::SppRef, T)
    #check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
            writedlm(sim, hcat("Running MORTALITY: ADULTS"))
    end

    old = findall(x -> (x.stage == "s" && x.age >= x.bankduration), plants)
    deleteat!(plants, old)

    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        writedlm(sim, hcat("Seeds dead on seed-bank:", length(old)))
    end

    dying = filter(x -> (x.stage == "s" && (rem(t,52) > x.seedoff || x.age > x.seedoff)), plants)
    
    for sp in unique(getfield.(dying, :sp))

    	Bm = seed_mfactor*sppref.b0mort[sp]*(sppref.seedmass[sp]^(-1/4))*exp(-aE/(Boltz*T))

	morts = factorized_seedproc("mortality", dying, Bm, sp, plants)
	deleteat!(plants, morts)
    end
    
end

function survive!(plants::Array{submodels.Plant,1}, t::Int, cK::Float64, K::Float64, settings::Dict{String, Any}, sppref::SppRef, landavail::BitArray{2},T, biomass_production::Float64, dying_stage::String)

    if dying_stage == "a"
        open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
            writedlm(sim, hcat("Running MORTALITY: ADULTS"))
        end
    elseif dying_stage == "j"
        open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
            writedlm(sim, hcat("Running MORTALITY: JUVENILES"))
        end
    end

    # Density-independent mortality
    # ------------------------------
    deaths = Int64[]
    mprob = 0
    Bm = 0
    b0mort = 0

    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        writedlm(sim, hcat("# total: ", length(plants),
		           "# seeds:", length(findall(x -> x.stage == "s", plants)),
		           "# juveniles:", length(findall(x -> x.stage == "j", plants)),
		           "# adults:", length(findall(x -> x.stage == "a", plants)),
		           "Above-ground vegetative weighing:", sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants), 0.00001))))
    end

    # old ones die
    old = findall( x -> (x.stage == dying_stage && x.age >= x.span), plants)
    deleteat!(plants, old)

    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        writedlm(sim, hcat("Dying of age:", length(old)))
    end

    # the rest of the individuals have a metabolic probability of dying. Seeds that are still in the mother plant cant die. If their release season is over, it is certain that they are not anymore, even if they have not germinated

    dying = findall(x -> x.stage == dying_stage, plants) # mortality function is run twice, focusing on juveniles or adults; it can only run once each, so seeds go with adults
    
    for d in dying
        # stages have different mortality factors
        if plants[d].stage == "j"
	    m_stage = juv_mfactor
        elseif plants[d].stage == "a"
	    m_stage = adult_mfactor
        else
	    error("Error with plant's stage assignment")
        end

        current_vegmass = plants[d].mass["leaves"] + plants[d].mass["stem"] + plants[d].mass["root"]
		Bm = plants[d].b0mort*m_stage*(current_vegmass^(-1/4))*exp(-aE/(Boltz*T))
        mprob = 1 - exp(-Bm)

        # unit test
        if mprob < 0
	    error("mprob < 0")
	    #mprob = 0
        elseif mprob > 1
	    error("mprob > 1")
        end

        if 1 == rand(Distributions.Bernoulli(mprob))
	    push!(deaths, d)
	    # check-points of life history processes and metabolic rates
	    open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
	        writedlm(sim, hcat(t, "death", plants[d].stage, plants[d].age))
	    end
	    open(abspath(joinpath(settings["outputat"],settings["simID"],"metaboliclog.txt")),"a") do sim
	        writedlm(sim, hcat(plants[d].stage, plants[d].age, Bm, mprob, "death"))
	    end
        end
    end
    deleteat!(plants, deaths) #delete the ones that are already dying due to mortality rate, so that they won't cramp up density-dependent mortality
    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        println(sim, "Density-independent mortality: ", length(deaths))
    end

    # Density-dependent mortality
    # ---------------------------
    deaths = Int64[] # reset before calculating density-dependent mortality
    Bm = 0
    mprob = 0
    b0mort = 0

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	println(sim, "Production before density-dependent mortality: $(sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants), 0.00001)))g; K = $K")
    end

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	writedlm(sim, hcat("# seeds:", length(findall(x -> x.stage == "s", plants)),
		           "# juveniles:", length(findall(x -> x.stage == "j", plants)),
		           "# adults:", length(findall(x -> x.stage == "a", plants)),
		           " Above-ground vegetative weighing:", sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants), 0.00001))))
    end

    # biomass of both juveniles and adults is used as criteria for production > K
    # but only only stage dies at each timestep
    biomass_plants = filter(x -> x.stage in ["j" "a"], plants) 

        # get coordinates of all occupied cells
	locs = getfield.(biomass_plants, :location)
	fullcells_indxs = findall(nonunique(DataFrame(hcat(locs))))

        if length(fullcells_indxs) > 0

            # check-point
            open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	        println(sim, "Number of shared cells: $(length(unique(locs[fullcells_indxs])))")
            end

	    for loc in unique(locs[fullcells_indxs])

		#find plants that are in the same grid
		plants_samecell = filter(x -> x.location == loc, biomass_plants)

		while sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants_samecell),0.00001)) > cK

		# check-point
            	open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	            println(sim, "Biomass in samecell $loc  with $(length(plants_samecell)) inds BEFORE dens-dep mort: 
		                  $(sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants_samecell),0.00001)))")
            	end
		
		    # get fitness of all species that are in the cell
                    sppgrid_fitness = Dict(sp => sppref.fitness[sp] for sp in unique(getfield.(plants_samecell, :sp)))

		    # check-point
            	    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	                 println(sim, "Sps sharing cell: $(keys(sppgrid_fitness)).")
            	    end
		    
		    # the species in the dictionnaries are over their respective cell
		    # carrying capacity, i.e., individuals must die. Therefore,
		    # loop through all species in the dictionary, killing accordingly.

		    for sp in keys(sppgrid_fitness)

                       	cK_sp = cK * (sppgrid_fitness[sp]/sum(collect(values(sppgrid_fitness))))

			# get individuals of sp in the current cell
                        plantssp_samecell = filter(x -> x.sp == sp, plants_samecell)
			    
                        # unit test
                        if (length(filter(x -> x.stage == "j", plantssp_samecell)) == 0 &&
                            length(filter(x -> x.stage == "a", plantssp_samecell)) == 0)
			    error("Cell carrying capacity overboard, but no juveniles or adults of $sp were detected")
		        end
			if(length(filter(x -> x.stage == "s", plantssp_samecell)) > 0)
			    error("Seeds being detected for density-dependent mortality")
		 	end

                        dying_plants = filter(x -> x.stage == dying_stage, plantssp_samecell)

                      	while (sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plantssp_samecell), 0.00001)) > cK_sp && length(dying_plants) > 0)

			    # check-point
            		    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	            	        println(sim, "Killing $sp inside grid-cell.")
            		    end
			    
                            # loop through smaller individuals (size instead of age, to keep things at a metabolic base)
                            masses = map(x -> (x.mass["leaves"]+x.mass["stem"]+x.mass["root"]), dying_plants)
		            dying = filter(x -> (x.mass["leaves"]+x.mass["stem"]+x.mass["root"]) == minimum(masses),
			    	    	   dying_plants)[1] #only one can be tracked down and killed at a time
					   		    #it is not possible to order `plants` by any field value

                            o = findall(x -> x.id == dying.id, plants)[1] #selecting "first" element changes the format into Int64, instead of native Array format returned by findall()

			    # check-point
 			    open(abspath(joinpath(settings["outputat"],settings["simID"],"eventslog.txt")),"a") do sim
		                writedlm(sim, hcat(t, "death-K-gridcell", plants[o].stage, plants[o].age))
                    	    end

			    deleteat!(plants, o)

			    # update control of while-loop
                            o_cell = findall(x -> x.id == dying.id, plantssp_samecell)[1] #selecting "first" element changes the format into Int64, instead of native Array format returned by findall()
                    	    deleteat!(plantssp_samecell, o_cell)
                            dying_plants = filter(x -> x.stage == dying_stage, plantssp_samecell)

                        end
                    end

		    # update of control of while-loop
		    biomass_plants = filter(x -> x.stage in ("j", "a"), plants)
       		    plants_samecell = filter(x -> x.location == loc, biomass_plants)
		    # check-point
            	    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	                 println(sim, "Biomass in samecell $loc with $(length(plants_samecell)) inds AFTER dens-dep mort: 
		                  $(sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants_samecell),0.00001)))")
                    end
                end
            end
    end

# check-point
open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
    println(sim, "Above-ground biomass after density-dependent mortality: $(sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants), 0.00001)))g; K = $K")
end

## Surviving ones get older: refilter, because some plants got deleted
if dying_stage == "a"
        surviving = findall(x -> ((x.stage == "s" && (rem(t,52) > x.seedoff || x.age > x.seedoff)) || x.stage == dying_stage), plants)
    else
        surviving = findall(x -> x.stage == dying_stage, plants) # mortality function is run twice, focusing on juveniles or adults; it can only run once each, so seeds go with adults
    end
for s in surviving
    plants[s].age += 1
end

end

"""
    shedflower!(plants, sppref, t, settings)
Plants loose their reproductive biomasses at the end of the reproductive season
 and 50% of biomass during winter.

"""
function shedflower!(plants::Array{submodels.Plant,1}, sppref::SppRef, t::Int, settings::Dict{String,Any})
    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        writedlm(sim, hcat("SHEDDING reproductive biomass/ WInter-die back"))
    end

    flowering = findall(x -> (x.mass["repr"] > 0 && rem(t,52) > x.floroff), plants)
    
    for f in flowering
	plants[f].mass["repr"] = 0.0
    end
    
end

"""
    winter_dieback!(plants, t)
In the last week of the year, all adult `plants` loose all of the biomass allocated to `leaves`
and reproductive (`repr`) structures. Biomass allocated to `stem` is decreased by a half, if not 
already at that value.
"""
function winter_dieback!(plants::Array{submodels.Plant,1}, t::Int)
    adults = findall(x -> (x.stage == "a"), plants)

    for a in adults
        plants[a].mass["leaves"] = 0.0
	plants[a].mass["repr"] = 0.0
	plants[a].mass["stem"] >= (0.5*plants[a].compartsize) ? plants[a].mass["stem"] = (0.5*plants[a].compartsize) : nothing
    end
end

"""
    destroyorgs!(plants)
Kill organisms that where in the lost habitat cells.

"""
function destroyorgs!(plants::Array{submodels.Plant,1}, landavail::BitArray{2}, settings::Dict{String,Any})
    kills = []
    for o in 1:length(plants)
	if landavail[plants[o].location...] == false
	    push!(kills,o)
	end
    end
    if length(kills) > 0 # trying to delete at index 0 generates an error
	deleteat!(plants, kills)
    end
end

"""
    manage!()
Juvenile and adult that are big enough, i.e., above-ground compartments have more than 50% its maximum value, have these compartments reduced to 50% of their biomass.
Mowing happens at most once a year, between August and September.
"""
function manage!(plants::Array{submodels.Plant,1}, t::Int64, management_counter::Int64, settings::Dict{String,Any})

    if management_counter < 1 || 1 == rand(Distributions.Bernoulli(manage_prob))

        # check-point
        open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	    writedlm(sim, hcat("Above-ground biomass before mowing =", sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants), 0.00001))))
        end

        mowed = findall(x -> (x.stage in ["j" "a"] &&
                           (x.mass["leaves"] >= (x.compartsize)^(3/4) || x.mass["stem"] >= 0.5*x.compartsize)), plants)

        for m in mowed
	    plants[m].mass["leaves"] >= (0.5*plants[m].compartsize)^(3/4) ? plants[m].mass["leaves"] = (0.5*plants[m].compartsize) : nothing
	    plants[m].mass["stem"] >= 0.5*plants[m].compartsize ? plants[m].mass["stem"] = (0.5*plants[m].compartsize) : nothing
	    plants[m].mass["repr"] = 0
        end

        management_counter += 1

        # check-point
        open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	    writedlm(sim, hcat("Biomass after MOWING =", sum(vcat(map(x -> (x.mass["leaves"]+x.mass["stem"]), plants), 0.00001))))
        end
    end
    return management_counter
end

"""
    landscape_init()
Create the initial landscape structure.
"""
function landscape_init(landpars::LandPars)

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

function landscape_init(landpars::NeutralLandPars)

    # convert matrix to BitArray (smaller than Bool)
    landavail = Bool.(landpars.initialland)

    # create landscape with same dimensions
    landscape = fill(Dict{String, Float64}(),
		     size(landavail))

    return landscape, landavail

end

"""
    updateenv!(landscape,t)
Update temperature and precipitation values according to the weekly input data (weekly means and ).
"""
function updateenv!(t::Int64, landpars::LandPars)

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

function updateenv!(t::Int64, landpars::NeutralLandPars)

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
destroyarea!()
Destroy proportion of habitat area according to input file. Destruction is simulated by making affected cells unavailable for germination and killing organisms in them.
"""
function destroyarea!(landpars::LandPars, landavail::Array{Bool,2}, settings::Dict{String,Any}, t::Int64)

    if settings["landmode"] == "artif"
	# DESTROY HABITAT
	# index of the cells still availble:
	available = findall(x -> x == true, landavail)
	# number of cells to be destroyed:
	loss = landpars.disturbarea/landpars.initialarea # just a proportion, unit doesn't matter
	lostarea = round(Int,loss*length(available), RoundUp)

	# Unity test
	if lostarea > length(available)
	    error("Destroying more than the available area.")
	end

	# go through landscape indexes of the first n cells from the still available
	for cell in available[1:lostarea]
	    landavail[cell] = false # and destroy them
	end

	#unity test
	open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"a") do sim
	    println(sim, "Number of destroyed cells: $lostarea")
	end

    elseif settings["landmode"] == "real"
	# rebuild the landscape according to shape file
	landscape = Array{Dict{String,Float64}, 2}

	for frag in collect(1:landpars.nfrags)

	    fragment = fill(Dict{String,Float64}(), landpars.flength[frag],landpars.flength[frag])

	    if frag == 1
		landscape = fragment #when empty, landscape cant cat with frag
	    else
		landscape = cat(3,landscape, frag)
	    end

	end

	landavail = fill(true,size(landscape))

    end
end

function destroyarea!(landpars::NeutralLandPars, landavail::BitArray{2}, settings::Dict{String,Any}, t::Int64)

    if settings["landmode"] == "artif"
	# DESTROY HABITAT
	# index of the cells still availble:
	available = findall(x -> x == true, landavail)
	# number of cells to be destroyed:
	loss = landpars.disturbland[:proportion][findall(landpars.disturbland[:td]==[t])[1]]
	lostarea = round(Int,loss*length(available), RoundUp)

	# Unity test
	if lostarea > length(available)
	    error("Destroying more than the available area.")
	end

	# go through landscape indexes of the first n cells from the still available
	for cell in available[1:lostarea]
	    landavail[cell] = false # and destroy them
	end

	#unity test
	open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"a") do sim
	    println(sim, "Number of destroyed cells: $lostarea")
	end

    elseif settings["landmode"] == "real"
	#TODO probably unnecessary
    end

    landscape = fill(Dict{String, Float64}(),
		     size(landavail))

    return landscape, landavail
end

"""
    fragment!()
This function is only called for simulating the fragmentation of an originally single landscape.
"""
function fragment!(landscape::Array{Dict{String,Float64},N} where N, settings::Dict{String,Any}, landpars::LandPars, plants::Array{submodels.Plant,1})

    newlandscape = []

    # Built fragmented landscape
    for frag in collect(1:landpars.nfrags)

	fragment = fill(Dict{String,Float64}(), landpars.flength[frag],landpars.flength[frag])

	if frag == 1
	    newlandscape = fragment #when empty, landscape cant cat with frag
	else
	    newlandscape = cat(3, newlandscape, fragment)
	end
    end

    # Resettle the individuals in the new landscape:
    # list all indexes of old landscape
    arrayidx = hcat(collect(ind2sub(landscape, findall(x -> x == x, landscape)))...)
    landscapeidx = [Tuple(arrayidx[x,:]) for x in 1:size(arrayidx, 1)]
    # create array to store the idx of the landscape cells that were not destroyed and sample old cells indexes to fill the new landscape
    remaincells = sample(landscapeidx, length(newlandscape); replace = false) # the size
    # reshap the cells in the same dimensions and sizes as the newlandscape, so the new indexes are correct
    idxholder = reshape(remaincells, size(newlandscape))
    # kill the organisms that have not remained in the new landscape configuration
    filter!(x -> x.location in idxholder, plants)
    # update their .location field
    newloc = []
    for o in plants
	newloc <- idxholder[findall(x->x == collect(plants[o].location), idxholder)]
	plants[o].location = newloc
    end

    landscape = newlandscape
    landavail = fill(true,size(landscape))

    return landscape, landavail

end

# method for 'continuous' landscape structure (not 3D)
function fragment!(landscape::Array{Dict{String, Float64}, N} where N, landavail::BitArray{2}, landpars::NeutralLandPars, t::Int64, tdist::Any)

    # convert matrix to BitArray (smaller than Bool)
    landavail = Bool.(landpars.disturbland)
    # create landscape with same dimensions
    landscape = fill(Dict{String, Float64}(),
		     size(landavail))

    return landscape, landavail
end

"""
    disconnect!()
Decreases the connectivity of an already fragmented landscape.
"""
function disconnect!()
end

end
