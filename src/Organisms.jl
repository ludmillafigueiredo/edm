"""
This module contains the data structures and functions to simulate plants, pollination regimes, and life cycle processes:
- allocate!(): biomass growth and resource allocation
- develop!(): maturation of juveniles that have reached their respective age of first flowering
- mate!(): pollination
- mkoffsrping!(): seed/clone production
- release!(): find seeds that can be dispersed
- disperse!(): set new locations for seeds being dispersed
- establish!(): see if seeds manage to germinate in current location
- shedd!(): decrease of reproductive biomass at the end of reproductive season, or winter die-back
- manage!(): decrease of vegetative and reproductive biomass due to mowing

Plants have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the range of search).
"""

module Organisms

using Distributions
using DataFrames
using JuliaDB
using DataValues
using StatsBase
using Fileprep

export SppRef, TraitRanges, Plant, initorgs, develop!, allocate!, mate!, mkoffspring!, microevolution!, disperse!, germinate, establish!, survive!, shedd!, manage!, destroyorgs!, release!

# Set up model constants
const Boltz = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
const aE = 0.63 #0.65 - eV Brown & Sibly MTE book chap 2
const µ_short = 1
const λ_short = 0.2
const µ_medium = 0.2
const λ_medium = 3
const µ_long = 1000
const λ_long = 100
const Q = 5

# upload glocal parameters
include("globalparameters_EDM.jl")

# Initial trait values is read from an input file and stored for reference in `sppref::SppRef`.
mutable struct SppRef
    sp_id::Array{String, 1}
    abund::Dict{String,Int}
    kernel::Dict{String,String}
    clonality::Dict{String,Bool}
    seedmass::Dict{String,Float64}
    maxmass::Dict{String,Float64}
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

# Minimal and maximal trait values, which control microevolution, are stored in `traitranges::TraitRanges`
mutable struct TraitRanges
    seedmass::Dict{String,Array{Float64,1}}
    maxmass::Dict{String,Array{Float64,1}}
    span::Dict{String,Array{Int64,1}}
    firstflower::Dict{String,Array{Int64,1}}
    floron::Dict{String,Array{Int64,1}}
    floroff::Dict{String,Array{Int64,1}}
    seednumber::Dict{String,Array{Int64,1}}
    seedon::Dict{String,Array{Int64,1}}
    seedoff::Dict{String,Array{Int64,1}}
    bankduration::Dict{String,Array{Int64,1}}
end

# Data structure holding individual trait values
mutable struct Plant
    id::String
    stage::String #e,j,a
    location::Tuple # (x,y)
    sp::String #sp id, easier to read
    kernel::String
    clonality::Bool
    #### Evolvable traits ####
    seedmass::Float64
    maxmass::Float64
    span::Int64
    firstflower::Int64
    floron::Int64
    floroff::Int64
    seednumber::Int64
    seedon::Int64
    seedoff::Int64
    bankduration::Int64
    b0grow::Float64
    b0germ::Float64
    b0mort::Float64
    #### State variables ####
    age::Int64 # control death when older than max. lifespan
    mass::Dict{String, Float64}
    mated::Bool
end

"""
initorgs(landavail, sppref,id_counter)

Initializes the organisms characterized in the input info stored in `sppref` and distributes them in the available landscape `landavail`. Stores theindividuals in the `orgs` array, which holds all organisms being simulated at any given time.

"""
function initorgs(landavail::BitArray{N} where N, sppref::SppRef, id_counter::Int, settings::Dict{String, Any}, K::Float64)

    plants = Plant[]

    for s in sppref.sp_id # all fragments are populated from the same species pool

        sp_abund = Int(round(((sppref.fitness[s]/sum(collect(values(sppref.fitness))))*K)/mean([sppref.seedmass[s],sppref.maxmass[s]]), RoundUp))
	# check-point
	open(abspath(joinpath(settings["outputat"], string(settings["simID"], "initialabundances.txt"))),"a") do sim
	    println(sim, "Initial abundance of $s: $sp_abund")
        end

	# create random locations
	XYs = hcat(rand(1:size(landavail,1), sp_abund),
		   rand(1:size(landavail,2), sp_abund))

	for i in 1:sp_abund

	    id_counter += 1 # update individual counter
	    minvalue = 1e-7 # Distribution.Uniform requires max > min

	    newplant = Plant(hex(id_counter),
			      rand(["a" "j" "s"]),
			      (XYs[i,1],XYs[i,2]),
			      s,
			      sppref.kernel[s],
			      sppref.clonality[s],
			      sppref.seedmass[s],
			      sppref.maxmass[s], #maxmass
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
			      Dict("veg" => 0.0, "repr" => 0.0),
			      false)
            
            ## initial biomass
	    if newplant.stage == "s"
		newplant.mass["veg"] = newplant.seedmass
		newplant.age = newplant.seedon + 1
	    elseif newplant.stage == "j"
		newplant.mass["veg"] = newplant.seedmass
		newplant.age = newplant.seedon + 4
	    elseif newplant.stage in ["a"]
		newplant.mass["veg"] = newplant.maxmass * 0.75
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
Calculates biomass gain according to the metabolic theory (`aE`, `Boltz` and `T` are necessary then). According to the week being simulated, `t` and the current state of the individual growing ( the biomass gained is

"""
function allocate!(plants::Array{Plant,1}, t::Int64, aE::Float64, Boltz::Float64, settings::Dict{String, Any},sppref::SppRef, T::Float64, biomass_production::Float64, K::Float64, growing_stage::String)
    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        writedlm(sim, hcat("Growth of", growing_stage))
    end

    growing = find(x->x.stage == growing_stage,plants)

    for o in growing

        b0grow = plants[o].b0grow
	
	B_grow = (b0grow*(plants[o].mass["veg"])^(-1/4))*exp(-aE/(Boltz*T))                                                           # only vegetative biomass fuels growth

	if plants[o].stage == "j"                                        # juveniles and small adults grow vegetative biomass only
            
	    new_mass = B_grow*(plants[o].maxmass - plants[o].mass["veg"])
	    plants[o].mass["veg"] += new_mass

        elseif (plants[o].stage == "a" &&
	        (plants[o].floron <= rem(t,52) < plants[o].floroff) &&
		plants[o].mass["veg"] >= 0.5*plants[o].maxmass)                                              # adults in their reprod. season invest in reproduction

            new_mass = B_grow*(plants[o].mass["veg"])

	    if haskey(plants[o].mass,"repr")
		plants[o].mass["repr"] += new_mass #sowingmass
	    else
		plants[o].mass["repr"] = new_mass #sowingmass
	    end
	elseif (plants[o].stage == "a" && 0.5*plants[o].maxmass < plants[o].mass["veg"] < plants[o].maxmass)
	    new_mass = B_grow*(plants[o].maxmass - plants[o].mass["veg"])
	    plants[o].mass["veg"] += new_mass

        end
	
    end

    # unity test
    masserror = find(x -> sum(collect(values(x.mass))) <= 0, plants)
    if length(masserror) > 0
        println("Plant: plants[masserror[1]]")
	error("Zero or negative values of biomass detected.")
    end

end

"""
    develop!()
Controls individual juvenile maturation.

"""
function develop!(plants::Array{Organisms.Plant,1}, sppref::SppRef, settings::Dict{String, Any}, t::Int)
    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        writedlm(sim, hcat("Maturation..."))
    end

    juvs = find(x->x.stage == "j",plants)

    for j in juvs
	if  plants[j].age >= plants[j].firstflower
	    plants[j].stage = "a"

	    # test
	    open(abspath(joinpath(settings["outputat"],settings["simID"],"eventslog.txt")),"a") do sim
		writedlm(sim, hcat(t, "maturation", plants[j].stage, plants[j].age))
	    end
	end
    end

end

"""
    mate!()
Calculate proportion of insects that reproduced (encounter?) and mark that proportion of the population with the `mated` label.
- visited: reduction in pollination service

"""
function mate!(plants::Array{Organisms.Plant,1}, t::Int, settings::Dict{String, Any}, scen::String, tdist::Any, remaining)

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        writedlm(sim, hcat("Pollination ..."))
    end

    ready = find(x-> x.stage == "a" && x.mass["repr"] > x.seedmass, plants)
    pollinated = []
    npoll = 0

    if length(ready) > 0 # check if there is anyone flowering

        # Scenarios were pollination is not species-specific
        # --------------------------------------------------

        if scen in ["indep" "equal"]
            
	    occupiedflwrs = rand(Distributions.Uniform(1e-4, 1e-2),1)[1]                                 # base number of pollinated flowes
	    defaultnpoll = rand(Distributions.Binomial(Int(ceil(length(ready) * occupiedflwrs)),0.5))[1]

	    # Determine number of individuals that get pollinated (species is not relevant)
	    if scen == "indep"
		npoll = defaultnpoll # Fishman & Hadany's proportion of visited flowers
		println("Scenario of INDEP pollination loss")

	    elseif scen == "equal" #all species lose pollination randomly (not species-specific)
		println("Scenario of EQUAL pollination loss")
		
		if t in tdist                                                          # calculate the amount of loss for the specified times
		    npoll = Int(ceil(defaultnpoll * remaining[find(tdist == [t])[1]]))
		else
		    npoll = defaultnpoll
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
After mating happened (marked in `reped`), calculate the amount of offspring

"""
function mkoffspring!(plants::Array{Organisms.Plant,1}, t::Int64, settings::Dict{String, Any},sppref::SppRef, id_counter::Int, landavail::BitArray{2}, T::Float64, traitranges::Organisms.TraitRanges)

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	writedlm(sim, hcat("Total number of individuals before REPRODUCTION:", length(plants)))
    end

    offspring = Plant[]
    non0sd = 1e-7
    seed_counter = 0

    # Separate sexually and asexually reproducing (mated status can change during the simulation and this would generate more clones)
    ferts = filter(x -> x.mated == true, plants)
    asexuals = filter(x -> x.mated == false && x.clonality == true && x.mass["repr"] > x.seedmass, plants)

    # Sexuallly produced offspring
    # ----------------------------
    for sp in unique(getfield.(ferts, :sp))

	sowing = find(x -> x.sp == sp || x.id in unique(getfield.(ferts, :id)), plants)

	spoffspringcounter = 0 #offspring is not output in the same file as adults and juveniles

	for s in sowing

	    seedmass = plants[s].seedmass
	    offs = div(0.5*plants[s].mass["repr"], seedmass)

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
		    for f in fieldnames(seed)
			if typeof(getfield(seed, f)) in [Int64 Float64]
			    if getfield(seed, f) < 0 || getfield(seed, f) == Inf
				error(f, " has value: ", getfield(seed, f), "Ind: ", seed.id, "sp: ", seed.sp)
			    end
			end
		    end

                    # State variables
		    seed.id = hex(id_counter)
		    seed.mass = Dict("veg" => seed.seedmass,
				       "repr" => 0.0)
                    seed.stage = "s"
                    seed.age = 0
                    seed.mated = false

                    # Trait microevolution
                    # --------------------
		    #seed.seedmass += rand(Distributions.Normal(0, abs(plants[s].seedmass-conspp.seedmass+non0sd)/6))[1]
                    seed.maxmass = (seed.maxmass+conspp.maxmass)/2+ rand(Distributions.Normal(0, abs(plants[s].maxmass-conspp.maxmass+non0sd)/6))[1]
		    seed.span = Int(round((seed.span+conspp.span)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].span-conspp.span+non0sd)/6))[1], RoundUp))
		    seed.firstflower = Int(round((seed.firstflower+conspp.firstflower)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].firstflower-conspp.firstflower+non0sd)/6))[1], RoundUp))
		    seed.floron = Int(round((seed.floron+conspp.floron)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].floron-conspp.floron+non0sd)/6))[1],RoundUp))
		    seed.floroff = Int(round((seed.floroff+conspp.floroff)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].floroff-conspp.floroff+non0sd)/6))[1],RoundUp))
		    seed.seednumber = Int(round((seed.number+conspp.seednumber)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].seednumber-conspp.seednumber+non0sd)/6))[1], RoundUp))
		    seed.seedon = Int(round((seed.seedon+conspp.seedon)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].seedon-conspp.seedon+non0sd)/6))[1],RoundUp))
		    seed.seedoff = Int(round((seed.seedoff+conspp.seedoff)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].seedoff-conspp.seedoff+non0sd)/6))[1],RoundUp))
		    seed.bankduration = Int(round((seed.bankduration+conspp.bankduration)/2, RoundUp)) + Int(round(rand(Distributions.Normal(0, abs(plants[s].bankduration-conspp.bankduration+non0sd)/6))[1],RoundUp))

		    # constrain microevolution: avoid trait changes that generate negative values (and also values that irrealistically high)

                    #if (seed.seedmass < traitranges.seedmass[seed.sp][1] || seed.seedmass > traitranges.seedmass[seed.sp][end])
		    #   seed.seedmass < traitranges.seedmass[seed.sp][1] ? seed.seedmass = traitranges.seedmass[seed.sp][1] :
		    #   seed.seedmass = traitranges.seedmass[seed.sp][end]
		    #end

		    if (seed.maxmass < traitranges.maxmass[seed.sp][1] || seed.maxmass > traitranges.maxmass[seed.sp][end])
			seed.maxmass < traitranges.maxmass[seed.sp][1] ? seed.maxmass = traitranges.maxmass[seed.sp][1] :
			    seed.maxmass = traitranges.maxmass[seed.sp][end]
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

# output seeds per species (file is initialized in main.jl)
open(abspath(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv")),"a") do seedfile
    writedlm(seedfile, hcat(t, sp, "s", "sex", spoffspringcounter))
end
end

# check-point
open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
    writedlm(sim, hcat("Total number of individuals after SEX:", length(plants)))
end


# Asexually produced offspring
# ----------------------------
for sp in unique(getfield.(asexuals, :sp))

    cloning = find(x -> x.sp == sp && x.id in unique(getfield.(asexuals, :id)) , plants) # find mothers cloning

    # start production counting
    spclonescounter = 0

    for c in cloning
	offs = div(0.5*plants[c].mass["repr"], plants[c].seedmass)

	# unity test
	if  plants[c].mass["repr"] <= 0
	    error("Negative reproductive biomass") #because `offs` is an integer, reproductive biomass should not become negative
	end

	if offs <= 0
	    continue
	else
	    # limit offspring production to the maximal number of seeds the species can produce
	    offs > plants[c].seednumber ? offs = plants[c].seednumber : offs

	    # update reproductive mass
	    plants[c].mass["repr"] -= (offs * plants[c].seedmass)

	    # get a copy of the mother, which the clones will look like
	    clonetemplate = deepcopy(plants[c])
	    clonetemplate.stage = "j" #clones have already germinated
	    clonetemplate.mass["veg"] = plants[c].maxmass*0.1
	    clonetemplate.mass["repr"] = 0.0

	    for o in offs

		clone = deepcopy(clonetemplate)

                # check if the new location is actually available before creating the clone
		# with cells at 1m2, clones stay in the same cell as the mother
                if checkbounds(Bool, landavail, clone.location[1], clone.location[2])

                    id_counter += 1
		    clone.id = hex(id_counter)

		    push!(plants, clone)
		    spclonescounter += 1
		end
	    end
	end
    end

    # output clones
    open(abspath(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv")),"a") do seedfile
	writedlm(seedfile, hcat(t, sp, "j", "asex", spclonescounter))
    end
end

# check-point
open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
    writedlm(sim, hcat("Total number of individuals after ASEX:", length(plants)))
end


return id_counter

end


"""
    release!()
Individuals being released in any given week are: in seed stage (=seed= & in their seed release period (seedon <= t <= seedoff for the species)

"""

function release!(plants::Array{Organisms.Plant,1}, t::Int, settings::Dict{String, Any},sppref::SppRef)

    seedsi = find(x -> x.stage == "s" && x.age == 0 && x.seedon <= rem(t,52) < x.seedoff, plants)

end

"""
    disperse!(landscape, plants, t, seetings, sppref,)
Seeds are dispersed.

"""
function disperse!(landavail::BitArray{2}, seedsi, plants::Array{Organisms.Plant, 1}, t::Int, settings::Dict{String, Any}, sppref::SppRef, landpars::Any, tdist::Any)#Setworld.LandPars)}

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        writedlm(sim, hcat("Dispersing ..."))
    end

    lost = Int64[]
    justdispersed = 0

    for d in seedsi                    # only seeds that have been released can disperse
	if plants[d].kernel == "short"
	    µ, λ = [µ_short λ_short]
	elseif plants[d].kernel == "medium"
	    µ, λ = [µ_medium λ_medium]
	elseif plants[d].kernel == "long"
	    µ, λ = [µ_long λ_long]
	elseif plants[d].kernel in ["medium-short", "short-medium"]
	    µ, λ = rand([[µ_short λ_short],
			 [µ_medium λ_medium]])
	elseif plants[d].kernel in ["medium-long", "long-medium"]
	    µ, λ = rand([[µ_long λ_long],
			 [µ_medium λ_medium]])
	elseif plants[d].kernel in ["long-short", "short-long"]
	    µ, λ = rand([[µ_short λ_short],
			 [µ_long λ_long]])
	elseif plants[d].kernel == "any"
	    µ,λ = rand([[µ_short λ_short],
			[µ_medium λ_medium],
			[µ_long λ_long]])
	else
	    error("Check dispersal kernel input for species $(plants[d].sp).")
	end

	dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ,λ),1)[1])

	# Find the cell to which it is dispersing
	θ = rand(Distributions.Uniform(0,2),1)[1]*pi
	xdest = plants[d].location[1] + dist*round(Int64, cos(θ), RoundNearestTiesAway)
	ydest = plants[d].location[2] + dist*round(Int64, sin(θ), RoundNearestTiesAway)

	if checkbounds(Bool, landavail, xdest, ydest) && landavail[xdest, ydest] == true  # Check if individual fall inside the habitat area, otherwise, discard it already
                                                                                          # checking the suitability first would make more sense but cant be done if cell is out of bounds
	    plants[d].location = (xdest,ydest)
	    justdispersed += 1

	else # if the new location is in an unavailable habitat or outside the landscape, the seed dies

	    push!(lost,d)
	    # check-point of life-history processes
	    open(abspath(joinpath(settings["outputat"],settings["simID"],"eventslog.txt")),"a") do sim
		writedlm(sim, hcat(t, "lost in dispersal", plants[d].stage, plants[d].age))
	    end
	end
    end
    
    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	println(sim, "Number of dispersing: $justdispersed")
    end
    
    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	writedlm(sim, hcat("Lost in dispersal:", length(lost)))
    end
    
    deleteat!(plants,lost)

end

"""
    establish!
Seed that have already been released (in the current time step, or previously - this is why `seedsi` does not limit who get to establish) and did not die during dispersal can establish.# only after release seed can establish. Part of the establishment actually accounts for the seed falling in an available cell. This is done in the dispersal() function, to avoid computing this function for individuals that should die anyway. When they land in such place, they have a chance of germinating (become seedlings - `j` - simulated by `germinate!`). Seeds that don't germinate stay in the seedbank, while the ones that are older than one year are eliminated.

"""
function establish!(plants::Array{Organisms.Plant,1}, t::Int, settings::Dict{String, Any}, sppref::SppRef, T::Float64, biomass_production::Float64, K::Float64)
    #REFERENCE: May et al. 2009
    
    establishing = find(x -> x.stage == "s", plants)

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	writedlm(sim, hcat("Seeds trying to ESTABLISH:", length(establishing)))
    end

    lost = Int64[]
    Bm = 0
    gprob = 0

    for o in establishing
	
        Bg = plants[o].b0germ*(plants[o].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
	gprob = 1-exp(-Bg)
	# test
	open(abspath(joinpath(settings["outputat"],settings["simID"],"metaboliclog.txt")),"a") do sim
	    writedlm(sim, hcat(plants[o].stage, plants[o].age, Bg, gprob, "germination"))
	end

	if gprob < 0
	    error("gprob < 0")
	elseif gprob > 1
	    error("gprob > 1")
	end

	germ = false
	if 1 == rand(Distributions.Bernoulli(gprob))
	    germ = true
	end

	if germ == true
	    plants[o].stage = "j"
	    plants[o].mass["veg"] = plants[o].seedmass

	    open(abspath(joinpath(settings["outputat"],settings["simID"],"eventslog.txt")),"a") do sim
		writedlm(sim, hcat(t, "germination", plants[o].stage, plants[o].age))
	    end

	end
    end

end

"""
    survive!(orgs, landscape)
Organism survival depends on total biomass, according to MTE rate. However, the proportionality constants (b_0) used depend on the cause of mortality: competition-related, where
plants in nogrwth are subjected to two probability rates

"""
function survive!(plants::Array{Organisms.Plant,1}, t::Int, cK::Float64, K::Float64, settings::Dict{String, Any}, sppref::SppRef, landavail::BitArray{2},T, biomass_production::Float64, dying_stage::String)

    # check-point
    if dying_stage == "a"
        open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
            writedlm(sim, hcat("Running MORTALITY: ADULTS & SEEDS"))
        end
    elseif dying_stage == "j"
        open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
            writedlm(sim, hcat("Running MORTALITY: JUVENILES"))
        end
    end

    # Density-independent mortality
    # ------------------------------
    deaths = Int64[]
    mprob = 0
    Bm = 0
    b0mort = 0

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        writedlm(sim, hcat("# total: ", length(plants),
		           "# seeds:", length(find(x -> x.stage == "s", plants)),
		           "# juveniles:", length(find(x -> x.stage == "j", plants)),
		           "# adults:", length(find(x -> x.stage == "a", plants)),
		           "vegetative weighing:", sum(vcat(map(x -> x.mass["veg"], plants), 0.00001))))
    end

    # old ones die
    old = find( x -> ((x.stage == "a" && x.age >= x.span)), plants) #|| (x.stage == "s" && x.age >= x.bankduration)), plants)
    deleteat!(plants, old)

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        writedlm(sim, hcat("Dying of age:", length(old)))
    end

    # the rest of the individuals have a metabolic probability of dying. Seeds that are still in the mother plant cant die. If their release season is over, it is certain that they are not anymore, even if they have not germinated
    if dying_stage == "a"
        dying = find(x -> ((x.stage == "s" && (rem(t,52) > x.seedoff || x.age > x.seedoff)) || x.stage == dying_stage), plants)
    else
        dying = find(x -> x.stage == dying_stage, plants) # mortality function is run twice, focusing on juveniles or adults; it can only run once each, so seeds go with adults
    end

    for d in dying

        # stages have different mortality factors
        if plants[d].stage == "s"
	    m_stage = seed_mfactor
        elseif plants[d].stage == "j"
	    m_stage = juv_mfactor
        elseif plants[d].stage == "a"
	    m_stage = adult_mfactor
        else
	    error("Error with plant's stage assignment") 
        end

        Bm = plants[d].b0mort*m_stage*(plants[d].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
        mprob = 1 - exp(-Bm)

        # unity test
        if mprob < 0
	    error("mprob < 0")
	    #mprob = 0
        elseif mprob > 1
	    error("mprob > 1")
        end

        if 1 == rand(Distributions.Bernoulli(mprob))
	    push!(deaths, d)
	    
	    # check-point of life history processes and metabolic rates
	    open(abspath(joinpath(settings["outputat"],settings["simID"],"eventslog.txt")),"a") do sim
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
	println(sim, "Production before density-dependent mortality: $(sum(vcat(map(x -> x.mass["veg"], plants), 0.00001)))g; K = $K")
    end

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	writedlm(sim, hcat("# seeds:", length(find(x -> x.stage == "s", plants)),
		           "# juveniles:", length(find(x -> x.stage == "j", plants)),
		           "# adults:", length(find(x -> x.stage == "a", plants)),
		           " vegetative weighing:", sum(vcat(map(x -> x.mass["veg"], plants), 0.00001))))
    end

    biomass_plants = filter(x -> x.stage in ["j" "a"], plants) #biomass of both juveniles and adults is used as criteria, but only only stage dies at each timestep

    if sum(vcat(map(x -> x.mass["veg"], biomass_plants), 0.00001)) > K

        # separate location coordinates and find all individuals that are in the same location as others (by compaing their locations with nonunique(only possible row-wise, not between tuples. This is the only way to get their indexes
	locs = map(x -> x.location, biomass_plants)
	fullcells_indxs = find(nonunique(DataFrame(hcat(locs))))

        # mortality intra grid cell first
	if length(fullcells_indxs) > 0

            # check-point
            open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	        println(sim, "number of occupied cells: $(length(unique(locs[fullcells_indxs])))")
            end

	    # loop through each grid cell, looking for the ones with biomass production over the cell carrying-capacity
	    for loc in unique(locs[fullcells_indxs])

		#find plants that are in the same grid
		samecell = filter(x -> x.location == loc, biomass_plants)

		while sum(vcat(map(x -> x.mass["veg"], samecell),0.00001)) > cK

		    # get fitness of all species that are in the cell
                    sppgrid_fitness = Dict(sp => sppref.fitness[sp] for sp in map(x -> x.sp, samecell))
                    # calculate the relative fitness, in relation to the other ones

		    # the species in the dictionnaries are over their respective cell carrying capacity, i.e., individuals must die. Therefore, loop through all species in the dictionary, killing accordingly.
		    for sp in keys(sppgrid_fitness)

                       	cK_sp = cK * (sppgrid_fitness[sp]/sum(collect(values(sppgrid_fitness))))
			# check-point

                        samecell_sp = filter(x -> x.sp == sp, samecell)

                        # unity test
                        if (length(filter(x -> x.stage == "j", samecell_sp)) == 0 &&
                            length(filter(x -> x.stage == "a", samecell_sp)) == 0)
			    error("Cell carrying capacity overboard, but no juveniles or adults of $sp were detected") 
		        end
			if(length(filter(x -> x.stage == "s", samecell_sp)) > 0)
			    error("Seeds being detected for density-dependent mortality")
		 	end	
                        
                        dying_plants = filter(x -> x.stage == dying_stage, samecell_sp)

                      	while (sum(vcat(map(x -> x.mass["veg"], samecell_sp), 0.00001)) > cK_sp && length(dying_plants) > 0)

                            # loop through smaller individuals (size instead of age, to keep things at a metabolic base)
                            masses = map(x -> x.mass["veg"], dying_plants)
		            dying = filter(x -> x.mass["veg"] == minimum(masses), dying_plants)[1] #but only one can be tracked down and killed at a time (not possible to order the `plants` array by any field value)
			    
                            o = find(x -> x.id == dying.id, plants)[1] #selecting "first" element changes the format into Int64, instead of native Array format returned by find()# check-point
 			    open(abspath(joinpath(settings["outputat"],settings["simID"],"eventslog.txt")),"a") do sim
		                writedlm(sim, hcat(t, "death-K-gridcell", plants[o].stage, plants[o].age))
                    	    end
		    	    deleteat!(plants, o)
          		    
			    # update control of while-loop  
                            o_cell = find(x -> x.id == dying.id, samecell_sp)[1] #selecting "first" element changes the format into Int64, instead of native Array format returned by find()
                    	    deleteat!(samecell_sp, o_cell)
                            dying_plants = filter(x -> x.stage == dying_stage, samecell_sp)
          		    
                        end
                    end

		    # update of control of while-loop
		    biomass_plants = filter(x -> x.stage in ("j", "a"), plants)
       		    samecell = filter(x -> x.location == loc, biomass_plants)
		    
                end
            end
	    
        else # unity test
            error("Biomass production over K in the landscape, but can't detect in which grid-cells this is happening.") # if production > K, it has to be detected at the grid-cell level.
        end
    end

# check-point
open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
    println(sim, "Production after density-dependent mortality: $(sum(vcat(map(x -> x.mass["veg"], plants), 0.00001)))g; K = $K")
end

## Surviving ones get older: some of them were not filtered above, so the ageing up need to be done separately to include all
for o in 1:length(plants)
    plants[o].age += 1
end

end

"""
    shedd!()
Plants loose their reproductive biomasses at the end of the reproductive season and 50% of biomass during winter.

"""
function shedd!(plants::Array{Organisms.Plant,1}, sppref::SppRef, t::Int, settings::Dict{String,Any})
    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        writedlm(sim, hcat("SHEDDING reproductive biomass/ WInter-die back"))
    end

    flowering = find(x -> (x.mass["repr"] > 0 && rem(t,52) > x.floroff), plants) #indexing a string returns a Char type, not String. Therefore, p must be Char ('').

    for f in flowering
	plants[f].mass["repr"] = 0
    end

    if (rem(t,52) == 51)

        adults = find(x -> (x.stage == "a"), plants)

	for a in adults
	    plants[a].mass["veg"] = 0.5*plants[a].mass["veg"]
	end
    end
end

"""
    destroyorgs!(plants)
Kill organisms that where in the lost habitat cells.

"""
function destroyorgs!(plants::Array{Organisms.Plant,1}, landavail::BitArray{2}, settings::Dict{String,Any})

    kills = []

    for o in 1:length(plants)

	if landavail[plants[o].location...] == false
	    push!(kills,o)
	end
    end

    if length(kills) > 0 # trying to delete at index 0 generates an error
	deleteat!(plants, kills)
    end
    # check-point
    #open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
    #println("Killed plants: $(length(kills))")
    #end
end

"""
    manage!()
Plants loose 20% of vegetative biomass and all of the reproductive biomass due to mowing. Mowing happens at most 3 times a year, between July and August.

"""
function manage!(plants::Array{Organisms.Plant,1}, t::Int64, management_counter::Int64, settings::Dict{String,Any})

    if management_counter < 1 || 1 == rand(Distributions.Bernoulli(0.5))

        # check-point
        open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	    writedlm(sim, hcat("Biomass before mowing =", sum(vcat(map(x -> x.mass["veg"], plants), 0.00001))))
        end
        
        adults = find(x -> (x.stage == "a"), plants)

        for a in adults
	    plants[a].mass["veg"] = 0.5*plants[a].mass["veg"]
	    plants[a].mass["repr"] = 0
        end

        management_counter += 1

        # check-point
        open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
	    writedlm(sim, hcat("Biomass after MOWING =", sum(vcat(map(x -> x.mass["veg"], plants), 0.00001))))
        end
    end
    return management_counter
end

"""
    pollination!()
Simulates plant-insect encounters and effective pollen transfer.

"""
#TODO find a not too cumbersome way of modelling pollen transfer: draw from Poisot's probability

# """
#     mate!()
# Insects reproduce if another one is found in the immediate vicinity.
# """
# function mate!(plant::Plantanism)
#     x, y, frag = plant.location #another plant of same sp should match the locations of focus
#     sp = plant.sp
#
#     # 1. check in the location field of plants array:
#     # 1.a inside same frag, look for locaions inside the squared area.
#     # 2. when matching, differentiate between autotrphsa and the rest
#     if plant.stage == "a"
#
#         for o in 1:length(plants) #look for partners
#             #TODO optimize indexation of field location in arrray
#             #TODO memory-wise, is it better to put all ifs together?
#             # 1:1 sex-ratio,
#             if frag == plants[o].location[3]
#                 if plants[o].location[1:2] in collect(Iterators.product(x-1:x+1,y-1:y+1))
#                     # check sp, self and already reproduced
#                     # if (sp == plants[o].sp && !(Base.isequal(plant, plants[o])) && plant.repr == false && plants[o].repr = false)
#                     #     #TODO add stochasticity
#                     #     plant.repr = true
#                     #     plants[o].repr = true
#                     #
#                     #     parents_genes = [plant.genotype, plants[o].repr]
#                     # end
#                 end
#             end
#         end
#
#     else
#         continue
#     end
#     return parents_genes #TODO check if it conflicts with modifying plants
# end

end