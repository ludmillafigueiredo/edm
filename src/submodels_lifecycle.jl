using Distributions
using DataFrames
using DataValues
using StatsBase
using DelimitedFiles

"""
    allocate!(orgs, t, setting, SPP_REF, T)
"""
function grow!(plants::Array{Plant,1}, t::Int64, settings::Dict{String, Any}, T::Float64, biomass_production::Float64, K::Float64, growing_stage::String)

    growing = filter(x-> x.stage == growing_stage, plants) 
    filter!(x-> x.stage != growing_stage, plants)

    flowering_ids = filter(x -> x.stage == "a" &&
    		    	        x.floron <= rem(t,52) < x.floroff &&
				(sum(values(x.mass))-x.mass["repr"]) >=
				0.5*(2*x.compartsize+x.compartsize^(3/4)),
			   growing) |>
	            x -> getfield.(x, :id)

    for sp in unique(getfield.(growing, :sp))
    	b0grow = SPP_REF.b0grow[sp]
    	growing_sp = filter(x->x.sp == sp, growing)
	map(x -> grow_allocate!(x, b0grow, flowering_ids), growing_sp)
	append!(plants, growing_sp)	
    end

    # unit test
    check_duplicates(plants)
    masserror = findall(x -> sum(collect(values(x.mass))) <= 0, plants)
    if length(masserror) > 0
        println("Plant: plants[masserror[1]]")
	error("Zero or negative values of biomass detected.")
    end
    
end

"""
    mature!(plants, settings, t)
Juveniles in `plants` older than their age of first flowering at time `t` become adults.
"""
function mature!(plants::Array{Plant,1}, settings::Dict{String, Any}, t::Int)

    nplants_before = length(plants)

    # find indexes of individuals that are ready to become adults
    juvs = filter(x-> (x.stage == "j" && x.age >= x.firstflower),plants)
    filter!(x -> !(x.id in getfield.(juvs, :id)), plants)
    # unit test
    if length(juvs) > 0
        if "a" in getfield.(juvs, :stage) || "s" in getfield.(juvs, :stage)
            error("Juvenile sorting for maturation is doing the opposite")
        end
    end
    setfield!.(juvs, :stage, "a")
    if "j" in getfield.(juvs, :stage)
            error("Maturation is not working")
        end
    append!(plants, juvs)
    if nplants_before != length(plants)
       error("Young adults got lost")
    end

    check_ages(plants)
    
    # unit test
    check_duplicates(plants)
    
end

"""
get_npoll
Calculate the number of individuals that get pollinated under each regime of pollinatoin.
Wind-pollinated species are not affected be pollination loss. Therefore, the number of pollinated individuals dependents only on the amount of flowering and on the default values of efficiency.
Scenarios of rdm and equal pollination loss have the same calculation, but for equal, the number is calculated per species, whereas for rdm, its the total number.
"""
function get_npoll(pollen_vector::String, poll_plants::Array{Plant,1}, poll_pars::PollPars, t::Int64)

    if pollen_vector == "wind"
        npoll = rand(Distributions.Binomial(Int(ceil(length(poll_plants)*WIND_DFLT)), WIND_EFFC))[1]
    else
	npoll_dflt = rand(Distributions.Binomial(Int(ceil(length(poll_plants)*VST_DFLT)),
						 INSCT_EFFC))[1]
	if poll_pars.scen == "indep"
     	    npoll = npoll_dflt 
	elseif poll_pars.scen in ["equal", "rdm"]
	    if t in poll_pars.regime.td
	        npoll = Int(ceil(npoll_dflt * poll_pars.regime[poll_pars.regime.td.== t,:remaining][1]))
	    else
		npoll = npoll_dflt
	    end
	elseif poll_pars.scen == "spec"
	    # pseudo
	end
    end

    return npoll

end

"""
pollinate!()
Mark plants and having pollinated (mated = true)
"""
function pollinate!(pollen_vector::String, flowering::Array{Plant,1}, plants::Array{Plant,1}, poll_pars::PollPars, t::Int64)

    poll_plants = filter(x -> occursin(pollen_vector, x.pollen_vector), flowering)
    npoll = get_npoll(pollen_vector, poll_plants, poll_pars, t)
    log_pollination(flowering, npoll, pollen_vector, t)
    pollinated = sample(poll_plants, npoll, replace = false, ordered = false)
    filter!(x -> !(x.id in getfield.(pollinated, :id)), flowering)
    setfield!.(pollinated, :mated, true)
    append!(plants, pollinated)
    
    # unit test
    check_duplicates(plants)
    
end

"""
    mate!(plants, t, settings)
Calculate proportion of `plants` that reproduced at time `t`, acording to pollination scenario `scen`, and mark that proportion of the population with the `mated` label.
"""
function mate!(plants::Array{Plant,1}, t::Int64, settings::Dict{String, Any}, poll_pars::PollPars)

    flowering = filter(x-> x.stage == "a" &&
    	    	       ALLOC_SEED*x.mass["repr"] > SPP_REF.seedmass[x.sp],
		       plants)
    if length(flowering) > 0 # check if there is anyone flowering

        # as with the other processes, it is easier to process plants separately from the main vector
        filter!(x -> !(x.id in getfield.(flowering, :id)), plants)
	
	pollinate!("wind", flowering, plants, poll_pars, t)
	
	if poll_pars.scen in ["equal", "rdm", "spec"] && t in poll_pars.regime.td
	    disturb_pollinate!(flowering, poll_pars, plants, t)
	else
	    pollinate!("insects", flowering, plants, poll_pars, t)
	end

	# re-insert the ones that did not get pollinated
	append!(plants, flowering)
    end

    # unit test
    check_duplicates(plants)
    
end

"""
microevolve!()
"""
function microevolve!(seed::Plant, fert::Plant, partner::Plant)

for trait in EVOLVABLE_TRAITS

    # Calculate microevolution
    # ------------------------
    new_traitvalue = mean([getfield(seed, trait), getfield(partner, trait)]) +
    		     rand(Normal(0, (abs(getfield(seed, trait)-getfield(partner, trait))+NOT_0)/6))[1]
    if trait != :compartsize
       new_traitvalue = Int(round(new_traitvalue, RoundUp))
    end
       
    # Constrain microevolution
    # ------------------------
    if new_traitvalue < getfield(TRAIT_RANGES, trait)[seed.sp][1]
       new_traitvalue = getfield(TRAIT_RANGES, trait)[seed.sp][1]
    elseif new_traitvalue > getfield(TRAIT_RANGES, trait)[seed.sp][end]
       new_traitvalue = getfield(TRAIT_RANGES, trait)[seed.sp][end]
    end   

    setfield!(seed, trait, new_traitvalue)

end
end

"""
    mkseeds!()
After mating happened (marked in `reped`), calculate the amount of offspring each individual produces, both sexually and assexually.
"""
function mkseeds!(plants::Array{Plant,1}, settings::Dict{String, Any}, T::Float64, t::Int64)

    # counters to keep track of offspring production
    seed_counter = 0

    # fertilized individuals
    ferts = filter(x -> x.mated == true, plants)

    # unit test
    if length(ferts) > length(findall(x -> x.stage == "a", plants))
       error("There are more fertilized individuals than adults")
    end
    
    for sp in unique(getfield.(ferts, :sp))

    	seedmass = SPP_REF.seedmass[sp]
	spoffspringcounter = 0 #offspring is not output in the same file as adults and juveniles
	
	ferts_sp = findall(x -> x.sp == sp && x.id in getfield.(ferts, :id), plants)
	
	for s in ferts_sp
	
	    offs = div(ALLOC_SEED*plants[s].mass["repr"], seedmass)
	    
	    if offs > 0
		
		# limit offspring production to the maximal number of seeds the species can produce
		offs > plants[s].seednumber ? offs = plants[s].seednumber : offs

                if t == 1 || rem(t,settings["output_freq"]) == 0
                    spoffspringcounter += offs
                    gather_offspring!(t, sp, "s", "sex", spoffspringcounter)
                end

		# Once it produces seeds, the plant looses the current "flowers" (repr. biomass)
		plants[s].mass["repr"] = 0

		# get a random parent
		partner = rand(ferts)

                for n in 1:offs

		    seed = deepcopy(plants[s])
		    seed_counter += 1
		    
		    # reassign state variables that dont evolve
		    seed.id = string(get_counter(), base = 16)
		    seed.mass = Dict("leaves" => 0.0,
		                     "stem" => 0.0,
				     "root" => seedmass,
				     "repr" => 0.0)
                    seed.stage = "s-in-flower"
                    seed.age = 0
                    seed.mated = false

		    if rand(Distributions.Binomial(1, 1 - plants[s].self_proba)) == 1
                        microevolve!(seed, plants[s], partner)
		    end
		    push!(plants, seed)

		end
	    end
        end

  end

    # unit test
    check_duplicates(plants)
    
end

"""
self_pollinate!()
"""
function self_pollinate!(plants::Array{Plant,1}, settings::Dict{String, Any}, t::Int64)

 selfers = filter(x-> x.stage == "a" &&
    	    	       ALLOC_SEED*x.mass["repr"] > 0.5*SPP_REF.seedmass[x.sp]x.seednumber &&
		       x.mated == false &&
		       x.self_failoutcross == true,
		       plants)|>
           x -> sample(x, Int(ceil(SELFING_PROBA*length(x))))

 for sp in unique(getfield.(selfers, :sp))

    	seedmass = SPP_REF.seedmass[sp]
	spoffspringcounter = 0 #offspring is not output in the same file as adults and juveniles
	
	selfers_sp = findall(x -> x.sp == sp && x.id in getfield.(selfers, :id), plants)
	
	for s in selfers_sp
	
	    offs = div(ALLOC_SEED*plants[s].mass["repr"], seedmass)
	    
	    if offs > 0
		
		# limit offspring production to the maximal number of seeds the species can produce
		offs > plants[s].seednumber ? offs = plants[s].seednumber : offs

                if t == 1 || rem(t,settings["output_freq"]) == 0
                    spoffspringcounter += offs
                    gather_offspring!(t, sp, "s", "sex", spoffspringcounter)
                end
                
		# Once it produces seeds, the plant looses the current "flowers" (repr. biomass)
		plants[s].mass["repr"] = 0

		for n in 1:offs

		    seed = deepcopy(plants[s])
		    
		    # reassign state variables that dont evolve
		    seed.id = string(get_counter(), base = 16)
		    seed.mass = Dict("leaves" => 0.0,
		                     "stem" => 0.0,
				     "root" => seedmass,
				     "repr" => 0.0)
                    seed.stage = "s-in-flower"
                    seed.age = 0
                    seed.mated = false

		    push!(plants, seed)

		end
	    end
        end
  end

  # unit test
  check_duplicates(plants)
  
end

"""
clone!()
Asexual reproduction.
"""
function clone!(plants::Array{Plant, 1}, settings::Dict{String, Any}, t::Int64)

    asexuals=filter(x -> x.mated==false && x.clonality==true && x.mass["repr"]>SPP_REF.seedmass[x.sp],
                    plants)

    for sp in unique(getfield.(asexuals, :sp))

        cloning = findall(x -> x.sp == sp && x.id in unique(getfield.(asexuals, :id)) , plants) # find mothers cloning

        # start species  production counting
        spclonescounter = 0

        for c in cloning

       	    if rand(Distributions.Binomial(1, 0.5)) == 1
                spclonescounter += 1
                clone = deepcopy(plants[c])
	        clone.stage = "j" #clones have already germinated
	        clone.age = 0
	        clone.mass["leaves"] = (plants[c].compartsize*0.1)^(3/4)
	        clone.mass["stem"] = plants[c].compartsize*0.1
	        clone.mass["root"] = plants[c].compartsize*0.1
	        clone.mass["repr"] = 0.0

	        clone.id = string(get_counter(), base=16)
	        push!(plants, clone)
                
       	    end
        end

        if t == 1 || rem(t,settings["output_freq"]) == 0
            gather_offspring!(t, sp, "s", "asex", spclonescounter)
        end
    end

    # unit test
    check_duplicates(plants)
    
end

"""
    disperse!(landscape, plants, t, seetings, SPP_REF,)
Seeds are dispersed.
`get_dest` is defined in `auxiliary.jl`
"""
function disperse!(landscape::BitArray{2},plants::Array{Plant, 1},t::Int,settings::Dict{String, Any},land_pars::Any)

    dispersing_idxs = findall(x -> (x.stage == "s-in-flower" && x.seedon <= rem(t,52) < x.seedoff),
                              plants)
    lost_idxs = Int64[]
    
    for i in dispersing_idxs
        
        kernel = split(plants[i].kernel, "-") |> collect |> rand
	
	# vectorized dispersal requires parameters to be vectorized too 
        dist = DISPERSAL_PARS[kernel].factor*
	rand(InverseGaussian(DISPERSAL_PARS[kernel].mu,
                             DISPERSAL_PARS[kernel].lambda))
	theta = rand(Uniform(0,2))*pi
	newloc = get_dest(plants[i].location, dist, theta)

        if newloc.suitable
            plants[i].stage = "s"
            plants[i].location = newloc.dest
        else
            push!(lost_idxs,i)
        end
    end

    deleteat!(plants, lost_idxs)
    if t == 1 || rem(t,settings["output_freq"]) == 0
        gather_event!(t, "lost in dispersal", "s", 0, length(lost_idxs))
        gather_event!(t, "dispersal", "s", 0, length(dispersing_idxs))
    end
    
end

"""
    establish!
Seed that have already been released (in the current time step, or previously - this is why `seedsi` does not limit who get to establish) and did not die during dispersal can establish in the grid-cell they are in. Germinated seeds mature to juveniles immediately. Seeds that don't germinate stay in the seedbank.
"""
function establish!(plants::Array{Plant,1}, t::Int, settings::Dict{String, Any},  T::Float64, biomass_production::Float64, K::Float64)

    # seeds that just dispersed try to establish and the ones that did not succeed previously
    # will get a chance again.
    establishing_idxs = findall(x -> x.stage == "s", plants)

    for i in establishing_idxs

    	if 1-exp(-(B0_GERM*(SPP_REF.seedmass[plants[i].sp]^(-1/4))*exp(-A_E/(BOLTZ*T)))) == 1
            plants[i].stage = "j"
            plants[i].age= 0
	end
    end

    if t == 1 || rem(t,settings["output_freq"]) == 0
        gather_event!(t, "germination", "j", 0, length(establishing_idxs))
    end
end

"""
    Plants density-independent mortality is calculated according to the metabolic theory.
Density-dependent mortality is calculated for cells with biomass over the grid-cell carrying capacity `C_K`. It kills smaller individuals until the biomass of the grid-cell is above `C_K` again.
Density-dependent mortality is not calculated for seeds. Their density-dependent mortality is calculated alongside adults just for convenience.
Both types of mortalities of adults and juveniles are calculated separately (as set by `dying_stage`).
    survive!(plants, t, settings, SPP_REF)
Calculate seed mortality, vectorized for all seeds of a same species.
"""
function die_seeds!(plants::Array{Plant,1}, settings::Dict{String, Any}, t::Int64, T::Float64)

    old = findall(x -> (x.stage == "s" && x.age > x.bankduration), plants)
    deleteat!(plants, old)
    
    dying = filter(x -> x.stage == "s", plants)
    n_deaths = 0

    n_plantsbefore = length(plants)
    
    for sp in unique(getfield.(dying, :sp))
    	dying_sp = filter(x -> x.sp == sp, plants)
    	Bm = SEED_MFACTOR*B0_MORT*(SPP_REF.seedmass[sp]^(-1/4))*exp(-A_E/(BOLTZ*T))
	death_idxs = vectorized_seedproc("mortality", dying_sp, Bm) |>
		    ids_deaths -> findall(x -> x.id in ids_deaths, plants)
	deleteat!(plants, death_idxs)
	n_deaths += length(death_idxs)
    end

    #unit test
    check_duplicates(plants)
    if n_plantsbefore - n_deaths != length(plants)
        error("Error in processing seed deaths")
    else
	log_sim("Seed mortality running smoooth")
    end

    if length(dying) > 0 && (t == 1 || rem(t,settings["output_freq"]) == 0)
        gather_event!(t, "death-indep", "s", mean(getfield.(dying, :age)), n_deaths)
    end
end

"""
die!()
"""
function die!(plants::Array{Plant, 1}, settings::Dict{String, Any}, T::Float64, dying_stage::String, t::Int64)

    old = findall( x -> (x.stage == dying_stage && x.age > x.span), plants)
    
    if length(old) > 0
       if dying_stage == "j"
           # unit test
	   log_age()
           error("Juvenile did not mature before reaching maximum span")
       else
	   deleteat!(plants, old)
       end 
    end

    # the rest of the individuals have a metabolic probability of dying.
    dying = filter(x -> x.stage == dying_stage, plants)
    filter!(x -> x.stage != dying_stage, plants)

    dead_ids, living_ids = survival(dying, T)
    dead_idxs = findall(x -> x.id in dead_ids, plants)
    deleteat!(plants, dead_idxs)
    living = filter(x -> x.id in living_ids, dying)
    append!(plants, living)

    # unit test
    check_duplicates(plants)

    if length(dying) > 0 && (t == 1 || rem(t,settings["output_freq"]) == 0)
        gather_event!(t, "death-indep", dying_stage, mean(getfield.(dying, :age)), length(dead_ids))
    end
end

function compete_die!(plants::Array{Plant,1}, t::Int, settings::Dict{String, Any},  landscape::BitArray{2}, T, dying_stage::String)

    
    # biomass of both juveniles and adults is used as criteria for check if  production > K
    # but only only stage dies at each timestep
    production_plants = filter(x -> x.stage in ["j" "a"], plants) 

    # get coordinates of all occupied cells
    locs = getfield.(production_plants, :location)
    fullcells_indxs = findall(nonunique(DataFrame(hcat(locs))))

    if length(fullcells_indxs) > 0

       for loc in unique(locs[fullcells_indxs])

           #find plants that are in the same grid
	   plants_cell = filter(x -> x.location == loc, production_plants)
  
           # get fitness of all species in the cell to simulate local competition
           sppcell_fitness = Dict(sp => SPP_REF.fitness[sp]
	                          for sp in unique(getfield.(plants_cell, :sp)))

           if sum(vcat(map(x->sum(values(x.mass))-x.mass["root"],plants_cell),NOT_0)) > C_K
	       for sp in keys(sppcell_fitness)
	           sort_die!(sp, sppcell_fitness, plants_cell, dying_stage, plants, settings, t)
	       end
           end
       end
    end

    # unit test
    check_duplicates(plants)
    
end

"""
    shedflower!(plants, SPP_REF, t, settings)
Plants loose their reproductive biomasses at the end of the reproductive season
 and 50% of biomass during winter.

"""
function shedflower!(plants::Array{Plant,1},  t::Int, settings::Dict{String,Any})

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
    manage!()
Juvenile and adult that are big enough, i.e., above-ground compartments have more than 50% its maximum value, have these compartments reduced to 50% of their biomass.
Mowing happens at most once a year, between August and September.
"""
function manage!(plants::Array{Plant,1}, t::Int64, management_counter::Int64, settings::Dict{String,Any})

    if management_counter < 1 || 1 == rand(Distributions.Bernoulli(MANAGE_PROB))

        mowed = findall(x -> (x.stage in ["j" "a"] &&
                           (x.mass["leaves"] >= (x.compartsize)^(3/4) || x.mass["stem"] >= 0.5*x.compartsize)), plants)

        for m in mowed
	    plants[m].mass["leaves"] >= (0.5*plants[m].compartsize)^(3/4) ? plants[m].mass["leaves"] = (0.5*plants[m].compartsize) : nothing
	    plants[m].mass["stem"] >= 0.5*plants[m].compartsize ? plants[m].mass["stem"] = (0.5*plants[m].compartsize) : nothing
	    plants[m].mass["repr"] = 0
        end

        management_counter += 1

    end
    
    return management_counter
    
end
