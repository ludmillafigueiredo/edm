using Distributions
using DataFrames
using DataValues
using StatsBase
using DelimitedFiles

"""
    allocate!(orgs, t, aE, Boltz, setting, SPP_REFERENCE, T)
"""
function allocate!(plants::Array{Plant,1}, t::Int64, aE::Float64, Boltz::Float64, settings::Dict{String, Any}, T::Float64, biomass_production::Float64, K::Float64, growing_stage::String)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Growth of $growing_stage")
    end

    growing = filter(x-> x.stage == growing_stage,plants)
	
    for sp in unique(getfield.(growing, :sp))
    
    	growing_sp = findall(x -> x.id in getfield.(growing, :id), plants)
    	b0grow = SPP_REFERENCE.b0grow[sp]
	
    for o in growing_sp

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
function develop!(plants::Array{Plant,1}, settings::Dict{String, Any}, t::Int)
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
function mate!(plants::Array{Plant,1}, t::Int, settings::Dict{String, Any}, scen::String, tdist::Any, remaining)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Pollination ...")
    end

    ready = findall(x-> x.stage == "a" && x.mass["repr"] > SPP_REFERENCE.seedmass[x.sp], plants)
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
    mkseeds!()
After mating happened (marked in `reped`), calculate the amount of offspring each individual produces, both sexually and assexually.
"""
function mkseeds!(plants::Array{Plant,1}, t::Int64, settings::Dict{String, Any}, id_counter::Int, landavail::BitArray{2}, T::Float64)

    # counters to keep track of offspring production
    seed_counter = 0

    # fertilized individuals
    ferts = filter(x -> x.mated == true, plants)

    # unit test
    if length(ferts) > length(findall(x -> x.stage == "a", plants))
       error("There are more fertilized individuals than adults")
    end
    
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Total number of individuals before REPRODUCTION: $length(plants)")
   	println(sim,
	"$(length(ferts)) fertilized in $(length(findall(x -> x.stage == "a", plants))) adults")
    end
    
    for sp in unique(getfield.(ferts, :sp))

    	seedmass = SPP_REFERENCE.seedmass[sp]
	spoffspringcounter = 0 #offspring is not output in the same file as adults and juveniles
	
	ferts_sp = findall(x -> x.sp == sp && x.id in getfield.(ferts, :id), plants)
	# check-point
    	open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	    println(sim, "Number of sowing: $(length(ferts_sp))")
        end

	for s in ferts_sp
	
	    offs = div(ALLOC_SEED*plants[s].mass["repr"], seedmass)
	    
	    if offs <= 0
		continue
	    else
		# limit offspring production to the maximal number of seeds the species can produce
		offs > plants[s].seednumber ? offs = plants[s].seednumber : offs

                spoffspringcounter += offs

		# Once it produces seeds, the plant looses the current "flowers" (repr. biomass)
		plants[s].mass["repr"] = 0

		# get a random parent
		partner = rand(ferts)

                for n in 1:offs

		    id_counter += 1

		    seed = deepcopy(plants[s])
		    seed_counter += 1

		    # unit test
		    for f in fieldnames(typeof(seed))
			if typeof(getfield(seed, f)) in [Int64 Float64]
			    if getfield(seed, f) < 0 || getfield(seed, f) == Inf
				error(f, " has value: ", getfield(seed, f), "Ind: ", seed.id, "sp: ", seed.sp)
			    end
			end
		    end

                    # reassign state variables that dont evolve
		    seed.id = string(id_counter, base = 16)
		    seed.mass = Dict("leaves" => 0.0,
		                     "stem" => 0.0,
				     "root" => seedmass,
				     "repr" => 0.0)
                    seed.stage = "s"
                    seed.age = 0
                    seed.mated = false

                    microevolve!(seed, plants[s], partner)

		    push!(plants, seed)

		end
	    end
        end

  	open(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv"),"a") do seedfile
    	    writedlm(seedfile, hcat(t, sp, "s", "sex", spoffspringcounter))
  	end
  end

    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        writedlm(sim, hcat("Total number of individuals after SEX:", length(plants)))
    end

    return id_counter
    
end

"""
microevolve!()
"""
function microevolve!(seed::Plant, fert::Plant, partner::Plant)

for trait in EVOLVABLE_TRAITS

    # Calculate microevolution
    # ------------------------
    new_traitvalue = mean([getfield(seed, trait), getfield(partner, trait)]) +
    		     rand(Normal(0, (abs(getfield(seed, trait)-getfield(partner, trait))+NOT_0)/6)[1]
    if trait == :compartsize
       new_traitvalue = Int(round(new_traitvalue, RoundUp))
    end
       
    # Constrain microevolution
    # ------------------------
    if new_traitvalue < getfield(TRAIT_RANGES, :trait)[seed.sp][1]
       new_traitvalue = getfield(TRAIT_RANGES, :trait)[seed.sp][1]
    end
    if new_traitvalue > getfield(TRAIT_RANGES, :trait)[seed.sp][end]
       new_traitvalue = getfield(TRAIT_RANGES, :trait)[seed.sp][end]
    end   

    setfield!(seed, trait, new_traitvalue)

end
end

"""
clone!()
Asexual reproduction.
"""
function clone!(plants)

    asexuals = filter(x -> x.mated == false && x.clonality == true && x.mass["repr"] > SPP_REFERENCE.seedmass[x.sp], plants)

    for sp in unique(getfield.(asexuals, :sp))

       cloning = findall(x -> x.sp == sp && x.id in unique(getfield.(asexuals, :id)) , plants) # find mothers cloning

       # start species  production counting
       spclonescounter = 0

       for c in cloning

       	if rand(Distributions.Binomial(1, 0.5)) == 1
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
    disperse!(landscape, plants, t, seetings, SPP_REFERENCE,)
Seeds are dispersed.

"""
function disperse!(landavail::BitArray{2}, plants::Array{Plant, 1}, t::Int, settings::Dict{String, Any},  landpars::Any, tdist::Any)

    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        writedlm(sim, hcat("Dispersing ..."))
    end

    lost = 0
    justdispersed = Plant[]

    #separating dispersing individuals facilitates indexing and assigning new field values
    dispersing = filter(x -> x.stage == "s" &&
    	       	 	x.age == 0 &&
			x.seedon <= rem(t,52) < x.seedoff, plants)
    filter!(x -> !(x.id in getfield.(dispersing, :id)), plants)
    
    for kernels in unique(getfield.(dispersing, :kernel))

    	dispersing_kernel = filter(x -> x.kernel == kernels, dispersing)
    	
	locs=[(idx=index,loc=dispersing_kernel[index].location) for index in 1:length(dispersing_kernel)]

        kernel = rand(split(kernels, "-") |> collect) #some have more than one kernel
        dists = dispersal_pars[kernel].factor*
	        rand(InverseGaussian(dispersal_pars[kernel].mu,dispersal_pars[kernel].lambda),length(locs))
	thetas = rand(Uniform(0,2), length(locs))*pi

	newlocs = map(get_dest, locs, dists, thetas)

	lost_idxs = filter(x->x.suitable==false, newlocs) |> x->getfield.(x,:idx)
	deleteat!(dispersing_kernel, lost_idxs)
	filter!(x->x.suitable==true, newlocs)
	
	setproperty!.(dispersing_kernel, :location, getfield.(newlocs, :dest))
	append!(justdispersed, dispersing_kernel)
	
	# check-point of life-history processes
	open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
	    for i in length(lost_idxs)
	        writedlm(sim, hcat(t, "lost in dispersal", "s", 0))
	    end
	end

	lost+=length(lost_idxs)
	
    end
    
    # update who actually dispersed
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	println(sim, "Number of successfull dispersals: $(length(justdispersed))")
    end
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	writedlm(sim, hcat("Lost in dispersal:", length(lost)))
    end
    
    return justdispersed
    
end

"""
    establish!
Seed that have already been released (in the current time step, or previously - this is why `seedsi` does not limit who get to establish) and did not die during dispersal can establish in the grid-cell they are in. Germinated seeds mature to juveniles immediately. Seeds that don't germinate stay in the seedbank.

"""
function establish!(justdispersed::Array{Plant,1}, plants::Array{Plant,1}, t::Int, settings::Dict{String, Any},  T::Float64, biomass_production::Float64, K::Float64)

    lost = Int64[]
    Bg = 0
    prob_g = 0

    establishing = filter(x -> x.stage == "s" && x.age >= 1, plants) |> x -> append!(justdispersed, x) #some young seeds might still be on mother plant
    filter!(x -> !(x.id in getfield.(establishing,:id)), plants)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	writedlm(sim, hcat("Seeds trying to ESTABLISH:", length(establishing)))
    end
    
    for sp in unique(getfield.(establishing, :sp))

    	Bg = SPP_REFERENCE.b0germ[sp]*(SPP_REFERENCE.seedmass[sp]^(-1/4))*exp(-aE/(Boltz*T))

	germinating_ids = factorized_seedproc("germination", establishing, Bg, sp, plants)

	germinating = filter(x -> x.id in germinating_ids, establishing)
	filter!(x -> !(x.id in germinating_ids), establishing)
	setproperty!.(germinating, :stage, "j")
	append!(plants, germinating)
    end
end

"""
    survive!(plants, t, cK, K, settings, SPP_REFERENCE, landavail, biomass_production, dying_stage)
Plants density-independent mortality is calculated according to the metabolic theory.
Density-dependent mortality is calculated for cells with biomass over the grid-cell carrying capacity `cK`. It kills smaller individuals until the biomass of the grid-cell is above `cK` again.
Density-dependent mortality is not calculated for seeds. Their density-dependent mortality is calculated alongside adults just for convenience.
Both types of mortalities of adults and juveniles are calculated separately (as set by `dying_stage`).
    survive!(plants, t, settings, SPP_REFERENCE)
Calculate seed mortality, vectorized for all seeds of a same species.
"""
function die_seeds!(plants::Array{Plant,1}, settings::Dict{String, Any}, t::Int64, T::Float64)
    old = findall(x -> (x.stage == "s" && x.age > x.bankduration), plants)
    deleteat!(plants, old)
    
    dying = filter(x -> x.stage == "s", plants)
    deaths = 0
    
    for sp in unique(getfield.(dying, :sp))
    	Bm = SEED_MFACTOR*B0_MORT*(SPP_REFERENCE.seedmass[sp]^(-1/4))*exp(-aE/(Boltz*T))
	death_idxs = vectorized_seedproc("mortality", dying, Bm, sp, plants) |>
		    ids_deaths -> findall(x -> x.id in ids_deaths, plants)
	deleteat!(plants, death_idxs)
	deaths += length(death_idxs)
    end

    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Seeds possibly dying: $(length(dying))")
	println(sim, "Dead on seed-bank: $(length(old))")
	println(sim, "Dead metabolic: $(length(deaths))")
    end

end

"""
die!()

"""
function die!(plants::Array{Plant, 1}, settings::Dict{String, Any}, T::Float64, dying_stage::String)

    old = findall( x -> (x.stage == dying_stage && x.age >= x.span), plants)
    deleteat!(plants, old)

    # the rest of the individuals have a metabolic probability of dying.
    dying = filter(x -> x.stage == dying_stage, plants)
    filter!(x -> x.stage != dying_stage, plants)

    dead_ids, living_ids = survival(dying, T)
    dead_idxs = findall(x -> x.id in dead_ids, plants)
    deleteat!(plants, dead_idxs)
    living = filter(x -> x.id in living_ids, dying)
    append!(plants, living)
    
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Density-independent mortality: $(length(dead_ids))")
        println(sim, "$(uppercase(dying_stage)) dying of age: $(length(old))")
    end

    
end

function compete_die!(plants::Array{Plant,1}, t::Int, cK::Float64, settings::Dict{String, Any},  landavail::BitArray{2}, T, dying_stage::String)

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        println(sim, "Above-ground biomass after density-dependent mortality:\n$(sum(vcat(map(x -> sum(values(x.mass))-x.mass["root"],
	filter(x -> x.stage in ["j", "a"], plants)), NOT_0)))g")
    end
    # biomass of both juveniles and adults is used as criteria for check if  production > K
    # but only only stage dies at each timestep
    production_plants = filter(x -> x.stage in ["j" "a"], plants) 

    # get coordinates of all occupied cells
    locs = getfield.(production_plants, :location)
    fullcells_indxs = findall(nonunique(DataFrame(hcat(locs))))

    if length(fullcells_indxs) > 0

       # check-point
       open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
           println(sim, "Number of shared cells: $(length(unique(locs[fullcells_indxs])))")
       end

	    for loc in unique(locs[fullcells_indxs])

		#find plants that are in the same grid
		plants_cell = filter(x -> x.location == loc, production_plants)

		# get fitness of all species in the cell to simulate local competition
                sppcell_fitness = Dict(sp => SPP_REFERENCE.fitness[sp]
		    		      	   for sp in unique(getfield.(plants_cell, :sp)))

		if sum(vcat(map(x->sum(values(x.mass))-x.mass["root"],plants_cell),NOT_0)) > cK

		    for sp in keys(sppcell_fitness)

                       	cK_sp = cK * (sppcell_fitness[sp]/sum(collect(values(sppcell_fitness))))

			prodplants_cell = filter(x -> x.stage in ["j", "a"] &&
					  	      x.sp == sp &&
						      x.location == loc,
					  	 plants)

			prodsp_cell = sum(vcat(map(x -> sum(values(x.mass))-x.mass["root"],prodplants_cell),NOT_0))
			
			if prodsp_cell > cK_sp

			   dying_plants = filter(x -> x.stage == dying_stage && x.sp == sp, plants_cell)
			   # order inds so smaller can be killed first with pop!()
                           dying_sorted = sort(dying_plants,
					       by = x -> sum(values(x.mass)), rev = true)

			   dying = Plant[]
			   while (prodsp_cell > cK_sp && length(dying_sorted) > 0)
			       # kill smallest
			       push!(dying, pop!(dying_sorted))

			       #update controls of while loop
			       prodplants_cell = filter(x -> x.stage in ["j", "a"] &&
			       		       	 	     x.sp == sp &&
							     x.location == loc,
					  	        plants)
			       prodsp_cell = sum(vcat(map(x -> sum(values(x.mass))-x.mass["root"],prodplants_cell),NOT_0))
                           end

			   dying_idxs = findall(x -> x.id in getfield.(dying, :id), plants)
			   deleteat!(plants, dying_idxs)

                        end

                    end
                end
            end
    end

    # check-point
    open(abspath(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt")),"a") do sim
        println(sim, "Above-ground biomass after density-dependent mortality:\n$(sum(vcat(map(x -> sum(values(x.mass))-x.mass["root"],
	filter(x -> x.stage in ["j", "a"], plants)), NOT_0)))g")
    end
    
end

"""
    shedflower!(plants, SPP_REFERENCE, t, settings)
Plants loose their reproductive biomasses at the end of the reproductive season
 and 50% of biomass during winter.

"""
function shedflower!(plants::Array{Plant,1},  t::Int, settings::Dict{String,Any})
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
function winter_dieback!(plants::Array{Plant,1}, t::Int)
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
function destroyorgs!(plants::Array{Plant,1}, landavail::BitArray{2}, settings::Dict{String,Any})
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
function manage!(plants::Array{Plant,1}, t::Int64, management_counter::Int64, settings::Dict{String,Any})

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

end

"""
    disturb!(landscape, landavail, plants, t, settings, landpars, tdist)
Change landscape structure and metrics according to the simulation scenario (`loss`, habitat loss, or `frag`, fragmentation)
"""
function disturb!(landscape::Array{Dict{String,Float64},N} where N, landavail::BitArray{2}, plants::Array{Plant,1}, t::Int64, settings::Dict{String,Any}, landpars::NeutralLandPars, tdist::Any)

    if settings["disturbtype"] == "loss"
         landscape, landavail = destroyarea!(landpars, landavail, settings, t)
    elseif settings["disturbtype"] == "frag"
         landscape, landavail = fragment!(landscape, landavail, landpars, t, tdist)
    end

    destroyorgs!(plants, landavail, settings)

    return landscape,landavail
end


"""
    fragment!()
This function is only called for simulating the fragmentation of an originally single landscape.
"""
function fragment!(landscape::Array{Dict{String,Float64},N} where N, settings::Dict{String,Any}, landpars::LandPars, plants::Array{Plant,1})

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

 end