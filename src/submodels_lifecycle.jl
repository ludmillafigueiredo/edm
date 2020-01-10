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
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Flowering: $(length(flowering_ids))")
    end

    log_abovegroundmass("before", "growth")
	
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
    # log process
    log_abovegroundmass("after", "growth")
    
end

"""
    mature!(plants, settings, t)
Juveniles in `plants` older than their age of first flowering at time `t` become adults.
"""
function mature!(plants::Array{Plant,1}, settings::Dict{String, Any}, t::Int)

    nplants_before = length(plants)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Maturation...")
    end

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
        npoll = rand(Distributions.Binomial(Int(ceil(length(poll_plants)*WIND_DFLT)),WIND_EFFC))[1]
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

# method for disturb_poll
function get_npoll(insects_plants::Array{Plant,1}, poll_pars::PollPars, t::Int64)

    npoll_dflt = rand(Distributions.Binomial(Int(ceil(length(insects_plants)*VST_DFLT)),
						 INSCT_EFFC))[1]
	if poll_pars.scen == "indep"
     	    npoll = npoll_dflt 
	elseif poll_pars.scen in ["equal", "rdm"]
	    if t in poll_pars.regime.td
	        npoll = Int(ceil(npoll_dflt * poll_pars.regime[poll_pars.regime.td.== t,
		                                               :remaining][1]))
	    else
		npoll = npoll_dflt
	    end
	elseif poll_pars.scen == "spec"
	    # pseudo
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
pollinate!()
Method used when pollination has been disturbed
"""
function pollinate!(insects_plants_sp::Array{Plant,1}, npoll::Int64, flowering::Array{Plant,1}, plants::Array{Plant,1}, poll_pars::PollPars, t::Int64)

    log_pollination(flowering, npoll, "insects-dist", t)
    pollinated = sample(insects_plants_sp, npoll, replace = false,
                        ordered = false)
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
    	    	       ALLOC_SEED*x.mass["repr"] > 0.5*SPP_REF.seedmass[x.sp]x.seednumber,
		       plants)
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Number of flowering plants: $(length(flowering)).")
    end

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
    
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Total number of individuals before REPRODUCTION: $(length(plants))")
   	println(sim,
	"$(length(ferts)) fertilized in $(length(findall(x -> x.stage == "a", plants))) adults")
    end
    
    for sp in unique(getfield.(ferts, :sp))

    	seedmass = SPP_REF.seedmass[sp]
	spoffspringcounter = 0 #offspring is not output in the same file as adults and juveniles
	
	ferts_sp = findall(x -> x.sp == sp && x.id in getfield.(ferts, :id), plants)
	# check-point
    	open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	    println(sim, "Number of $sp sowing: $(length(ferts_sp))")
        end

	for s in ferts_sp
	
	    offs = div(ALLOC_SEED*plants[s].mass["repr"], seedmass)
	    
	    if offs > 0
		
		# limit offspring production to the maximal number of seeds the species can produce
		offs > plants[s].seednumber ? offs = plants[s].seednumber : offs

                spoffspringcounter += offs
		open(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv"),"a") do seedfile
    	            writedlm(seedfile, hcat(t, sp, "s", "sex", spoffspringcounter))
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
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        writedlm(sim, hcat("Total number of individuals after SEX:", length(plants)))
    end    
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
	# check-point
    	open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	    println(sim, "Number of selfers: $(length(selfers_sp))")
        end

	for s in selfers_sp
	
	    offs = div(ALLOC_SEED*plants[s].mass["repr"], seedmass)
	    
	    if offs > 0
		
		# limit offspring production to the maximal number of seeds the species can produce
		offs > plants[s].seednumber ? offs = plants[s].seednumber : offs

                spoffspringcounter += offs
		open(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv"),"a") do seedfile
    	            writedlm(seedfile, hcat(t, sp, "s", "self", spoffspringcounter))
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
       	clone = deepcopy(plants[c])
	clone.stage = "j" #clones have already germinated
	clone.age = 0
	clone.mass["leaves"] = (plants[c].compartsize*0.1)^(3/4)
	clone.mass["stem"] = plants[c].compartsize*0.1
	clone.mass["root"] = plants[c].compartsize*0.1
	clone.mass["repr"] = 0.0

	clone.id = string(get_counter(), base=16)
	push!(plants, clone)
	spclonescounter += 1
            # check-point of life-history processes
	    open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
	        writedlm(sim, hcat(t, "clonal reproduction", "a", plants[c].age))
	    end
       	end
       end

       # output clones
        open(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv"),"a") do seedfile
	writedlm(seedfile, hcat(t, sp, "j", "asex", spclonescounter))
    	end
    end

    # unit test
    check_duplicates(plants)
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        writedlm(sim, hcat("Total number of individuals after ASEX:", length(plants)))
    end
    
end

"""
    disperse!(landscape, plants, t, seetings, SPP_REF,)
Seeds are dispersed.
`get_dest` is defined in `auxiliary.jl`
"""
function disperse!(landscape::BitArray{2},plants::Array{Plant, 1},t::Int,settings::Dict{String, Any},land_pars::Any)

    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        writedlm(sim, hcat("Dispersing ..."))
    end

    lost = 0
    justdispersed = Plant[]

    #separating dispersing individuals facilitates indexing and assigning new locations
    dispersing = filter(x -> (x.stage == "s-in-flower" &&
			x.seedon <= rem(t,52) < x.seedoff), plants)
    n_dispersing = length(dispersing)
    filter!(x -> !(x.id in getfield.(dispersing, :id)), plants)

    # check-point of life-history processes
    open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
        for i in 1:length(dispersing)
            writedlm(sim, hcat(t, "dispersal", "s", mean(getfield.(dispersing, :age))))
        end
    end
    
    for kernels in unique(getfield.(dispersing, :kernel))

    	dispersing_kernel = filter(x -> x.kernel == kernels, dispersing)

	# use indexes of individuals inside the vector with the ones being currently dispersed
	# facilitates keeping track of successes and losses in a vectorized manner.
	locs=[(idx=index,loc=dispersing_kernel[index].location)
	      for index in 1:length(dispersing_kernel)]

	# some dispersal sydromes are composed of several kernels.
	# one is randomly chosen for a particular instance of dispersal at time step.
        kernel = rand(split(kernels, "-") |> collect)
	
	# vectorized dispersal requires parameters to be vectorized too 
        dists = DISPERSAL_PARS[kernel].factor*
	        rand(InverseGaussian(DISPERSAL_PARS[kernel].mu,DISPERSAL_PARS[kernel].lambda),
		     length(locs))	     
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

    # unit test
    if lost +  length(justdispersed) != n_dispersing
       error("Number of successes and fails in dispersal do not match.")
    end
    check_duplicates(plants)
    
    # for now, just change the stage: these seeds are no longer in a flower
    # these individuals will be put back into the main vector `plants` after establishment (next)
    setproperty!.(justdispersed, :stage, "s")
 
    # check point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	println(sim, "Number of successfull dispersals: $(length(justdispersed))")
    	println(sim, "Lost in dispersal: $lost")
    end

    return justdispersed
    
end

"""
    establish!
Seed that have already been released (in the current time step, or previously - this is why `seedsi` does not limit who get to establish) and did not die during dispersal can establish in the grid-cell they are in. Germinated seeds mature to juveniles immediately. Seeds that don't germinate stay in the seedbank.
"""
function establish!(justdispersed::Array{Plant,1}, plants::Array{Plant,1}, t::Int, settings::Dict{String, Any},  T::Float64, biomass_production::Float64, K::Float64)

    # seeds that just dispersed try to establish and the ones that did not succeed previously
    # will get a chance again.
    establishing = filter(x -> x.stage == "s", plants) |> x -> append!(justdispersed, x)
    n_establishing = length(establishing)
    filter!(x -> !(x.id in getfield.(establishing,:id)), plants)
    
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	println(sim, "Seeds trying to establish: $n_establishing")
    end

    germinations = 0
    non_germinations = 0
    
    for sp in unique(getfield.(establishing, :sp))

    	establishing_sp = filter(x->x.sp==sp, establishing)
		
    	Bg = B0_GERM*(SPP_REF.seedmass[sp]^(-1/4))*exp(-A_E/(BOLTZ*T))
	
	germinated_ids = vectorized_seedproc("germination", establishing_sp, Bg)
	germinated = filter(x -> x.id in germinated_ids, establishing_sp)
        setproperty!.(germinated, :stage, "j")
	setproperty!.(germinated, :age, 0)
	append!(plants, germinated)
	germinations += length(germinated)

	# update the seeds that did not germinate and will go back into the main vector
	non_germinated = filter(x -> !(x.id in germinated_ids), establishing_sp)
    	append!(plants, non_germinated)
	non_germinations += length(non_germinated)

	# check-point of life-history processes
	open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
            for i in 1:length(germinated)
	        writedlm(sim, hcat(t, "germination", "j", mean(getfield.(germinated, :age))))
            end
	end
	
    end

    # check-point: can probably go, becuase the unit test is enough
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
	println(sim, "Seeds failing to establish: $non_germinations")
	println(sim, "Seeds germinating: $germinations")
    end

    # unit test
    check_duplicates(plants)
    if germinations + non_germinations != n_establishing
       error("Number of successes and fails in establishment do not match.")
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
    deaths = 0

    n_plantsbefore = length(plants)
    
    for sp in unique(getfield.(dying, :sp))
    	dying_sp = filter(x -> x.sp == "sp", plants)
    	Bm = SEED_MFACTOR*B0_MORT*(SPP_REF.seedmass[sp]^(-1/4))*exp(-A_E/(BOLTZ*T))
	death_idxs = vectorized_seedproc("mortality", dying_sp, Bm) |>
		    ids_deaths -> findall(x -> x.id in ids_deaths, plants)
	deleteat!(plants, death_idxs)
	deaths += length(death_idxs)
    end

    #unit test
    check_duplicates(plants)
    if n_plantsbefore - deaths != length(plants)
        error("Error in processing seed deaths")
    else
	log_sim("Seed mortality running smoooth")
    end

    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Seeds possibly dying: $(length(dying))")
	println(sim, "Dead on seed-bank: $(length(old))")
	println(sim, "Dead metabolic: $(length(deaths))")
    end
    open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
        for i in 1:deaths
            writedlm(sim, hcat(t, "death-indep", "s", mean(getfield.(dying, :age))))
        end
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
    
    # check-point
    open(joinpath(settings["outputat"],settings["simID"],"checkpoint.txt"),"a") do sim
        println(sim, "Density-independent mortality: $(length(dead_ids))")
        println(sim, "$(uppercase(dying_stage)) dying of age: $(length(old))")
    end
    open(joinpath(settings["outputat"],settings["simID"],"eventslog.txt"),"a") do sim
        for i in 1:length(dead_ids)
            writedlm(sim, hcat(t, "death-indep", dying_stage, mean(getfield.(dying, :age))))
        end
    end
    
end

function compete_die!(plants::Array{Plant,1}, t::Int, settings::Dict{String, Any},  landscape::BitArray{2}, T, dying_stage::String)

    # log process
    log_abovegroundmass("before", "density-depedent mortality")
    
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
    # log process
    log_abovegroundmass("after", "density-dependent mortality")
    
end

"""
    shedflower!(plants, SPP_REF, t, settings)
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
    manage!()
Juvenile and adult that are big enough, i.e., above-ground compartments have more than 50% its maximum value, have these compartments reduced to 50% of their biomass.
Mowing happens at most once a year, between August and September.
"""
function manage!(plants::Array{Plant,1}, t::Int64, management_counter::Int64, settings::Dict{String,Any})

    if management_counter < 1 || 1 == rand(Distributions.Bernoulli(MANAGE_PROB))

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
