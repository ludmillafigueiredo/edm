"""
                                    This module contains the

                                    Organisms have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the range of search).
                                    """
module Organisms

using Distributions
using DataFrames
using JuliaDB
using DataValues
using StatsBase
#using Setworld
using Fileprep

export OrgsRef_normal, OrgsRef_unif, TraitRanges, Organism, initorgs, develop!, allocate!, mate!, mkoffspring!, disperse!, germinate, establish!, survive!, shedd!, destroyorgs!, release!

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

# Initial organisms parametrization is read from an input file and stored in OrgsRef
mutable struct OrgsRef_normal
    sp_id::Array{String, 1}
    abund::Dict{String,Int}
    kernel::Dict{String,String}
    clonality::Dict{String,Bool}
    seedmass::Dict{String,Float64}
    maxmass::Dict{String,Float64}
    span_mean::Dict{String,Float64}
    span_sd::Dict{String,Float64}
    firstflower_mean::Dict{String,Float64}
    firstflower_sd::Dict{String,Float64}
    floron_mean::Dict{String,Float64}
    floron_sd::Dict{String,Float64}
    floroff_mean::Dict{String,Float64}
    floroff_sd::Dict{String,Float64}
    seednumber_mean::Dict{String,Float64}
    seednumber_sd::Dict{String,Float64}
    seedon_mean::Dict{String,Float64}
    seedon_sd::Dict{String,Float64}
    seedoff_mean::Dict{String,Float64}
    seedoff_sd::Dict{String,Float64}
    bankduration_mean::Dict{String,Float64}
    bankduration_sd::Dict{String,Float64}
    b0grow_mean::Dict{String,Float64}
    b0grow_sd::Dict{String,Float64}
    b0germ_mean::Dict{String,Float64}
    b0germ_sd::Dict{String,Float64}
    b0mort_mean::Dict{String,Float64}
    b0mort_sd::Dict{String,Float64}
end

mutable struct OrgsRef_unif
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
    floron_min::Dict{String,Float64}
    floron_max::Dict{String,Float64}
    floroff_min::Dict{String,Float64}
    floroff_max::Dict{String,Float64}
    seednumber_min::Dict{String,Float64}
    seednumber_max::Dict{String,Float64}
    seedon_min::Dict{String,Float64}
    seedon_max::Dict{String,Float64}
    seedoff_min::Dict{String,Float64}
    seedoff_max::Dict{String,Float64}
    bankduration_min::Dict{String,Float64}
    bankduration_max::Dict{String,Float64}
    b0grow_min::Dict{String,Float64}
    b0grow_max::Dict{String,Float64}
    b0germ_min::Dict{String,Float64}
    b0germ_max::Dict{String,Float64}
    b0mort_min::Dict{String,Float64}
    b0mort_max::Dict{String,Float64}
end

mutable struct TraitRanges
    seedmass::Dict{String,UnitRange{Int64}}
    maxmass::Dict{String,UnitRange{Int64}}
    span::Dict{String,UnitRange{Int64}}
    firstflower::Dict{String,UnitRange{Int64}}
    floron::Dict{String,UnitRange{Int64}}
    floroff::Dict{String,UnitRange{Int64}}
    seednumber::Dict{String,UnitRange{Int64}}
    seedon::Dict{String,UnitRange{Int64}}
    seedoff::Dict{String,UnitRange{Int64}}
    bankduration::Dict{String,UnitRange{Int64}}
    b0grow::Dict{String,UnitRange{Int64}}
    b0germ::Dict{String,UnitRange{Int64}}
    b0mort::Dict{String,UnitRange{Int64}}
end

mutable struct Organism
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
                        initorgs(landavail, orgsref,id_counter)

                        Initializes the organisms characterized in the input info stored in `orgsref` and distributes them in the available landscape `landavail`. Stores theindividuals in the `orgs` array, which holds all organisms being simulated at any given time.

                        """
function initorgs(landavail::BitArray{N} where N, orgsref, id_counter::Int, settings::Dict{String, Any})

    orgs = Organism[]
    
    for s in orgsref.sp_id # all fragments are populated from the same species pool

        # create random locations
        XYs = hcat(rand(1:size(landavail,1),orgsref.abund[s]),
                   rand(1:size(landavail,2),orgsref.abund[s]))

        for i in 1:orgsref.abund[s]

            id_counter += 1 # update individual counter
	    minvalue = 1e-7 # Distribution.Normal requires sd > 0, and Distribution.Uniform requires max > min

            if settings["traitdist"] == "uniform"
                neworg = Organism(hex(id_counter),
                                  rand(["a" "j" "e"]),
                                  (XYs[i,1],XYs[i,2]),
                                  s,
                                  orgsref.kernel[s],
                                  orgsref.clonality[s],
                                  orgsref.seedmass[s],
			          orgsref.maxmass[s], #maxmass
			          Int(round(rand(Distributions.Uniform(orgsref.span_min[s], orgsref.span_max[s] + minvalue),1)[1], RoundUp)),
			          Int(round(rand(Distributions.Uniform(orgsref.firstflower_min[s], orgsref.firstflower_max[s] + minvalue),1)[1], RoundUp)),
                                  Int(round(rand(Distributions.Uniform(orgsref.floron_min[s],orgsref.floron_max[s] + minvalue),1)[1], RoundUp)),
                                  Int(round(rand(Distributions.Uniform(orgsref.floroff_min[s],orgsref.floroff_max[s] + minvalue),1)[1], RoundUp)),
			          0, #seed number
			          Int(round(rand(Distributions.Uniform(orgsref.seedon_min[s],orgsref.seedon_max[s] + minvalue),1)[1], RoundUp)),
                                  Int(round(rand(Distributions.Uniform(orgsref.seedoff_min[s],orgsref.seedoff_max[s] + minvalue),1)[1], RoundUp)),
			          Int(round(rand(Distributions.Uniform(orgsref.bankduration_min[s],orgsref.bankduration_max[s] + minvalue),1)[1], RoundUp)),
			          rand(Distributions.Uniform(orgsref.b0grow_min[s],orgsref.b0grow_max[s] + minvalue),1)[1],
			          rand(Distributions.Uniform(orgsref.b0germ_min[s],orgsref.b0germ_max[s] + minvalue),1)[1],
			          rand(Distributions.Uniform(orgsref.b0mort_min[s],orgsref.b0mort_max[s] + minvalue),1)[1],
                                  0, #age
                                  Dict("veg" => 0.0, "repr" => 0.0), #mass
                                  false) #mated

                ## Set conditional traits and variables
                # weekly number of seeds
	        nseeds = rand(Distributions.Uniform(orgsref.seednumber_min[s],orgsref.seednumber_max[s] + minvalue),1)[1]
	        
            elseif settings["traitsdist"] == "normal"
                neworg = Organism(hex(id_counter),
                                  rand(["a" "j" "e"]),
                                  (XYs[i,1],XYs[i,2]),
                                  s,
                                  orgsref.kernel[s],
                                  orgsref.clonality[s],
                                  orgsref.seedmass[s],
			          orgsref.maxmass[s], #maxmass
			          Int(round(rand(Distributions.Normal(orgsref.span_mean[s], orgsref.span_sd[s] + minvalue),1)[1], RoundUp)),
			          Int(round(rand(Distributions.Normal(orgsref.firstflower_mean[s], orgsref.firstflower_sd[s] + minvalue),1)[1], RoundUp)),
                                  Int(round(rand(Distributions.Normal(orgsref.floron_mean[s],orgsref.floron_sd[s] + minvalue),1)[1], RoundUp)),
                                  Int(round(rand(Distributions.Normal(orgsref.floroff_mean[s],orgsref.floroff_sd[s] + minvalue),1)[1], RoundUp)),
			          0, #seed number
			          Int(round(rand(Distributions.Normal(orgsref.seedon_mean[s],orgsref.seedon_sd[s] + minvalue),1)[1], RoundUp)),
                                  Int(round(rand(Distributions.Normal(orgsref.seedoff_mean[s],orgsref.seedoff_sd[s] + minvalue),1)[1], RoundUp)),
			          Int(round(rand(Distributions.Normal(orgsref.bankduration_mean[s],orgsref.bankduration_sd[s] + minvalue),1)[1], RoundUp)),
			          rand(Distributions.Normal(orgsref.b0grow_mean[s],orgsref.b0grow_sd[s] + minvalue),1)[1],
			          rand(Distributions.Normal(orgsref.b0germ_mean[s],orgsref.b0germ_sd[s] + minvalue),1)[1],
			          rand(Distributions.Normal(orgsref.b0mort_mean[s],orgsref.b0mort_sd[s] + minvalue),1)[1],
                                  0, #age
                                  Dict("veg" => 0.0, "repr" => 0.0), #mass
                                  false) #mated

                ## Set conditional traits and variables
                # weekly number of seeds
	        nseeds = rand(Distributions.Normal(orgsref.seednumber_mean[s],orgsref.seednumber_sd[s] + minvalue),1)[1]
	        
            end

            # adjustments to traits than depend on other values
            if (neworg.floroff-neworg.floron+1) < 0
	        error("floroff - floron <0")
	    end
            neworg.seednumber = Int(round(nseeds/(neworg.floroff-neworg.floron+1), RoundUp))
            
            ## initial biomass
            if neworg.stage == "e"                  
                neworg.mass["veg"] = neworg.seedmass
                neworg.age = 1
            elseif neworg.stage in ["j"]
                neworg.mass["veg"] = neworg.maxmass * 0.25
                neworg.age = 26
	    elseif neworg.stage in ["a"]
                neworg.mass["veg"] = neworg.maxmass * 0.5
                neworg.age = 26
            else
                error("Check individual stages.")
            end

push!(orgs, neworg)

end

end

return orgs, id_counter

end

"""
                    allocate!(orgs, t, aE, Boltz, setting, orgsref, T)
                    Calculates biomass gain according to the metabolic theory (`aE`, `Boltz` and `T` are necessary then). According to the week being simulated, `t` and the current state of the individual growing ( the biomass gained is
                                    """
function allocate!(orgs::Array{Organism,1}, t::Int64, aE::Float64, Boltz::Float64, settings::Dict{String, Any},orgsref,T::Float64)
    #1. Initialize storage of those that dont growi and will have higher prob of dying (later)
    nogrowth = Int64[]

    growing = find(x->(x.stage in ["a" "j"]),orgs)
    for o in growing

        b0 = orgs[o].b0grow
        
        #only vegetative biomass helps growth
        B_grow = b0*(orgs[o].mass["veg"])^(-1/4)*exp(-aE/(Boltz*T))
        
        if orgs[o].stage == "j"
            # juveniles grow vegetative biomass only
	    new_mass = B_grow*(orgs[o].maxmass - orgs[o].mass["veg"])
            orgs[o].mass["veg"] += new_mass 
        elseif orgs[o].stage == "a" &&
            (orgs[o].floron <= rem(t,52) < orgs[o].floroff) &&
            (sum(collect(values(orgs[o].mass))) >= 0.5*(orgs[o].maxmass))
            # adults in their reproductive season and with enough weight, invest in reproduction
            #sowingmass = (5.5*(10.0^(-2)))*((orgs[o].mass["veg"]/(orgs[o].floroff-orgs[o].floron + 1) + new_mass)^0.95)
            new_mass = B_grow*(orgs[o].mass["veg"])
            if haskey(orgs[o].mass,"repr")
                orgs[o].mass["repr"] += new_mass #sowingmass 
            else
                orgs[o].mass["repr"] = new_mass #sowingmass
            end
        elseif orgs[o].stage == "a" && orgs[o].mass["veg"] < orgs[o].maxmass
            # adults that have not yet reached maximum size can still grow vegetative biomass, independently of the season
            new_mass = B_grow*(orgs[o].maxmass - orgs[o].mass["veg"])  
            orgs[o].mass["veg"] += new_mass 
        end            
    end
    return nogrowth
end

"""
                    develop!()
                    Controls individual juvenile maturation.
                    """
function develop!(orgs::Array{Organism,1}, orgsref)
    juvs = find(x->x.stage == "j",orgs)

    for j in juvs
        if  orgs[j].age >= orgs[j].firstflower
            # If an individual grows quite fast, it is more vigorous, and should transfer it to adult fecundity. The only variable capable of transfering this property is the weigh, which, combined with the MTE rate, makes it  generate more offspring
            orgs[j].stage = "a"
	end
    end
    
end

"""
                mate!()
                Calculate proportion of insects that reproduced (encounter?) and mark that proportion of the population with the `mated` label.
                - visited: reduction in pollination service
                """
function mate!(orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any}, scen::String, tdist::Any, remaining)

    ready = find(x-> x.stage == "a" && x.mass["repr"] > x.seedmass, orgs) # TODO find those with higher reproductive mas than the mean nb of seeds * seed mass.
    pollinated = []
    npoll = 0

    if length(ready) > 0 # check if there is anyone flowering
        # Scenarios were pollination is not species-specific
        if scen in ["indep" "equal"]  # calculation of number of pollinated individuals is different, but the actual pollination (non species-specific) is

            # Base number of pollinated flowes
	    occupiedflwrs = rand(Distributions.Uniform(1e-4, 1e-2),1)[1]
            defaultnpoll = rand(Distributions.Binomial(Int(ceil(length(ready) * occupiedflwrs)),0.5))[1] 
            
            # Determine number of individuals that get pollinated (species is not relevant)
            if scen == "indep"
                npoll = defaultnpoll # Fishman & Hadany's proportion of visited flowers
                println("Scenario of INDEP pollination loss")

            elseif scen == "equal" #all species lose pollination randomly (not species-specific)
                println("Scenario of EQUAL pollination loss")
                # calculate the amount of loss for the specified times
                if t in tdist
                    npoll = Int(ceil(defaultnpoll * remaining[find(tdist == [t])[1]]))
                else
                    npoll = defaultnpoll
                end
            end

            # Non species-specific pollination
            # check if any should actually be pollinated
            if npoll > 0
                # who are they
                pollinated = sample(ready, npoll, replace = false, ordered = true)
                # pollinate
                for p in pollinated
                    orgs[p].mated = true
                end
            elseif npoll < 0
                error("Negative number of plants being pollinated.")
            end

        elseif scen in ["rdm" "spec"] #not yet tested

            # Determine which species will loose pollination
            if scen == "rdm"
                # randomly pick plant species that will loose pollination (n = nspppoll) at a given timestep and find their number
                # pseudo:
                # rdmly pick species:
                # spppoll = unique(getfields(orgs, :sp)) |> sample(, nspppoll)
            elseif scen == "spec" #not yet tested
                # from a list of loss pollinators, find the plant species (in the interaction matrix) that will loose pollination at a given timestep
            end

            # Species-specific pollination
            # check if any should be pollinated
            if npoll == 0

            elseif npoll > 1
                # who are they
                pollinated = sample(ready, npoll, replace = false, ordered = true)
                # pollinate
                for p in pollinated
                    orgs[p].mated = true
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
function mkoffspring!(orgs::Array{Organisms.Organism,1}, t::Int64, settings::Dict{String, Any},orgsref, id_counter::Int, landavail::BitArray{2}, T::Float64, traitranges::Organisms.TraitRanges)

    offspring = Organism[]
    non0sd = 1e-7
    
    # Sexually produced offspring
    ferts = filter(x -> x.mated == true, orgs)

    for sp in unique(getfield.(ferts, :sp))

        sowing = find(x -> x.sp == sp || x.id in unique(getfield.(ferts, :id)), orgs)

        spoffspringcounter = 0

        for s in sowing

            seedmass = orgs[s].seedmass
            offs = div(0.5*orgs[s].mass["repr"], seedmass)

            if offs <= 0
                continue
            else
                # limit offspring production to the maximal number of seeds the species can produce
                offs > orgs[s].seednumber ? offs = orgs[s].seednumber : offs

                # update available biomass for reproduction
                orgs[s].mass["repr"] -= (offs * seedmass)

                # unity test
                if  orgs[s].mass["repr"] <= 0
                    error("Negative reproductive biomass")
                end

                # count species offspring for output
                spoffspringcounter += offs

                # get another parent
                conspp = orgs[rand(find(x -> x.sp == sp && x.stage == "a", orgs))] 

                for n in 1:offs

                    id_counter += 1 # update individual counter

                    embryo = deepcopy(orgs[s])

                    # unity test
                    for f in fieldnames(embryo)
		    	if typeof(getfield(embryo, f)) in [Int64 Float64]
		    	    if getfield(embryo, f) < 0 || getfield(embryo, f) == Inf
			        error(f, " has value: ", getfield(embryo, f), "Ind: ", embryo.id, "sp: ", embryo.sp)
			    end
			end
		    end

                    #
		    newvalue_seedmass = rand(Distributions.Normal(orgs[s].seedmass, abs(orgs[s].seedmass-conspp.seedmass+non0sd)/6))[1]
		    embryo.maxmass = rand(Distributions.Normal(orgs[s].maxmass, abs(orgs[s].maxmass-conspp.maxmass+non0sd)/6))[1]
		    embryo.span = Int(round(rand(Distributions.Normal(orgs[s].span, abs(orgs[s].span-conspp.span+non0sd)/6))[1], RoundUp))
                   embryo.firstflower = Int(round(rand(Distributions.Normal(orgs[s].firstflower, abs(orgs[s].firstflower-conspp.firstflower+non0sd)/6))[1], RoundUp))
		    embryo.floron = Int(round(rand(Distributions.Normal(orgs[s].floron, abs(orgs[s].floron-conspp.floron+non0sd)/6))[1],RoundUp))
		    embryo.floroff = Int(round(rand(Distributions.Normal(orgs[s].floroff, abs(orgs[s].floroff-conspp.floroff+non0sd)/6))[1],RoundUp))
		    embryo.seednumber = Int(round(rand(Distributions.Normal(orgs[s].seednumber, abs(orgs[s].seednumber-conspp.seednumber+non0sd)/6))[1], RoundUp))
		    embryo.seedon = Int(round(rand(Distributions.Normal(orgs[s].seedon, abs(orgs[s].seedon-conspp.seedon+non0sd)/6))[1],RoundUp))
		    embryo.seedoff = Int(round(rand(Distributions.Normal(orgs[s].seedoff, abs(orgs[s].seedoff-conspp.seedoff+non0sd)/6))[1],RoundUp))
		    embryo.bankduration = Int(round(rand(Distributions.Normal(orgs[s].bankduration, abs(orgs[s].bankduration-conspp.bankduration+non0sd)/6))[1], RoundUp))
		    embryo.b0grow = rand(Distributions.Normal(orgs[s].b0grow, abs(orgs[s].b0grow-conspp.b0grow+non0sd)/6))[1]
		    embryo.b0germ = rand(Distributions.Normal(orgs[s].b0germ, abs(orgs[s].b0germ-conspp.b0germ+non0sd)/6))[1]
		    embryo.b0mort = rand(Distributions.Normal(orgs[s].b0mort, abs(orgs[s].b0mort-conspp.b0mort+non0sd)/6))[1]
		    
                    # constrain values: avoid to trait changes that generates negative values (and also values that get too high)	   
                    
                    embryo.seedmass in traitranges.seedmass[embryo.sp] ? continue :
		   	embryo.seedmass < traitranges.seedmass[embryo.sp][1] ? embryo.seedmass == 0.95*orgsref.seedmass[embryo.sp] :
			embryo.seedmass == 1.05*orgsref.seedmass[embryo.sp]
		    embryo.maxmass in traitranges.maxmass[embryo.sp] ? continue :
		   	embryo.maxmass < traitranges.maxmass[embryo.sp][1] ? embryo.maxmass == 0.95*orgsref.seedmass[embryo.sp] :
			embryo.maxmass == 1.05*orgsref.seedmass[embryo.sp]
	            embryo.span in traitranges.span[embryo.sp] ? continue :
		   	embryo.span < traitranges.span[embryo.sp][1] ? embryo.span == 0.95*orgsref.span_min[embryo.sp] :
			embryo.span == 1.05*orgsref.span_max[embryo.sp]
	            embryo.firstflower in traitranges.firstflower[embryo.sp] ? continue :
		   	embryo.firstflower < traitranges.firstflower[embryo.sp][1] ? embryo.firstflower == 0.95*orgsref.firstflower_min[embryo.sp] :
			embryo.firstflower == 1.05*orgsref.firstflower_max[embryo.sp]
	            embryo.floron in traitranges.floron[embryo.sp] ? continue :
		   	embryo.floron < traitranges.floron[embryo.sp][1] ? embryo.floron == 0.95*orgsref.floron_min[embryo.sp] :
			embryo.floron == 1.05*orgsref.floron_max[embryo.sp]
	            embryo.floroff in traitranges.floroff[embryo.sp] ? continue :
		   	embryo.floroff < traitranges.floroff[embryo.sp][1] ? embryo.floroff == 0.95*orgsref.floroff_min[embryo.sp] :
		        embryo.floroff == 1.05*orgsref.floroff_max[embryo.sp]
                    embryo.seednumber in traitranges.seednumber[embryo.sp] ? continue :
		    	embryo.seednumber < traitranges.seednumber[embryo.sp][1] ? embryo.seednumber == 0.95*orgsref.seednumber_min[embryo.sp] :
			embryo.seednumber == 1.05*orgsref.seednumber_max[embryo.sp]
		    embryo.seedon in traitranges.seedon[embryo.sp] ? continue :
	   	        embryo.seedon < traitranges.seedon[embryo.sp][1] ? embryo.seedon == 0.95*orgsref.seedon_min[embryo.sp] :
			embryo.seedon == 1.05*orgsref.seedon_max[embryo.sp]
		    embryo.seedoff in traitranges.seedoff[embryo.sp] ? continue :
		   	embryo.seedoff < traitranges.seedoff[embryo.sp][1] ? embryo.seedoff == 0.95*orgsref.seedof_min[embryo.sp] :
			embryo.seedoff == 1.05*orgsref.seedoff_max[embryo.sp]
		    embryo.bankduration in traitranges.bankduration[embryo.sp] ? continue :
		   	embryo.bankduration < traitranges.bankduration[embryo.sp][1] ? embryo.bankduration == 0.95*orgsref.bankduration_min[embryo.sp] :
			embryo.bankduration == 1.05*orgsref.bankduration_max[embryo.sp]
		    embryo.b0grow in traitranges.b0grow[embryo.sp] ? continue :
		   	embryo.b0grow < traitranges.b0grow[embryo.sp][1] ? embryo.b0grow == 0.95*orgsref.b0grow_min[embryo.sp] :
			embryo.b0grow == 1.05*orgsref.b0grow_max[embryo.sp]
                    embryo.b0germ in traitranges.b0germ[embryo.sp] ? continue :
                        embryo.b0germ < traitranges.b0germ[embryo.sp][1] ? embryo.b0germ == 0.95*orgsref.b0germ_min[embryo.sp] :
                        embryo.b0germ == 1.05*orgsref.b0germ_max[embryo.sp]
                    embryo.b0mort in traitranges.b0mort[embryo.sp] ? continue :
                        embryo.b0mort < traitranges.b0mort[embryo.sp][1] ? embryo.b0mort == 0.95*orgsref.b0mort_min[embryo.sp] :
                        embryo.b0mort == 1.05*orgsref.b0mort_max[embryo.sp]

# set embryos state variables
embryo.id = hex(id_counter) 
embryo.mass = Dict("veg" => embryo.seedmass,
                   "repr" => 0.0)
embryo.stage = "e"
embryo.age = 0
embryo.mated = false

push!(offspring, embryo)
#unity test
#open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a")do sim
#println("Pushed new org ", embryo, " into offspring")
#end
end
orgs[s].mated = false # after producing seeds in a week, the plant will only do it again in the next week if it gets pollinated again
end 
end

# output seeds per species (file is initialized in main.jl)
open(abspath(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv")),"a") do seedfile
    writedlm(seedfile, hcat(t, sp, "e", "sex", spoffspringcounter))
end
end

# Asexually produced offspring
asexuals = filter(x -> x.mated == false && x.clonality == true && x.mass["repr"] > x.seedmass, orgs) #find which the individuals the can reproduce assexually and then go through them, by species

for sp in unique(getfield.(asexuals, :sp))
    cloning = find(x -> x.sp == sp && x.id in unique(getfield.(asexuals, :id)) , orgs)

    # start production counting
    spclonescounter = 0

    for c in cloning # mothers cloning
        offs = div(0.5*orgs[c].mass["repr"], orgs[c].seedmass)

        # unity test
        if  orgs[c].mass["repr"] <= 0
            error("Negative reproductive biomass") #because offs is an integer, reproductive biomass should not become negative
        end

        if offs <= 0
            continue
        else
            # limit offspring production to the maximal number of seeds the species can produce
            offs > orgs[c].seednumber ? offs = orgs[c].seednumber : offs

            # count clonal production  
            spclonescounter += offs

            # update reproductive mass
            orgs[c].mass["repr"] -= (offs * orgs[c].seedmass)

            # get a copy of the mother, which the clones will look like 
            clonetemplate = deepcopy(orgs[c])
            clonetemplate.stage = "j" #clones have already germinated
            clonetemplate.mass["veg"] = orgs[c].maxmass*0.1
            clonetemplate.mass["repr"] = 0.0

            for o in offs

                clone = deepcopy(clonetemplate)
                # check if the new location is actually available before creating the clone
                clone.location = (clonetemplate.location[1] + Int(round(rand(Distributions.Uniform(-1,1))[1])),
                                  clonetemplate.location[2] + Int(round(rand(Distributions.Uniform(-1,1))[1]))) # clones are spread in one of the neighboring cells - or in the same as the mother

                if checkbounds(Bool, landavail, clone.location[1], clone.location[2])
                    # actually start the new individual
                    id_counter += 1
                    clone.id = hex(id_counter)                  
                    
		    push!(offspring, clone)
		end
            end
        end
    end

    # output clones per species (file is initialized in main.jl)
    open(abspath(joinpath(settings["outputat"],settings["simID"],"offspringproduction.csv")),"a") do seedfile
        writedlm(seedfile, hcat(t, sp, "j", "asex", spclonescounter))
    end
end

append!(orgs, offspring)

#unity test
open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"a") do sim
    println(sim, "Offspring = " ,length(offspring), "and adults = ", length(find(x -> x.stage == "a", orgs)))
end

return id_counter

end

"""
                                    release!()
                                    Probably need to call disperse here, because not all "e"s in orgs are released at the same time.
                                    """
function release!(orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any},orgsref)

    # Individuals being released in any given week are: in embryo stage (=seed= & in their seed release period (seedon <= t <= seedoff for the species)
    seedsi = find(x -> x.stage == "e" && x.age == 0 && x.seedon <= rem(t,52) < x.seedoff, orgs)
    # using a condition "outside" orgs might not work. This condition with orgsref only works because orgsref always has the sp names of x as keys in the dictionnary. If presented with a key that it does do contain, it throws an error.

    return seedsi
end

"""
                    disperse!(landscape, orgs, t, seetings, orgsref,)
                    Seeds are dispersed.
                    """

function disperse!(landavail::BitArray{2}, seedsi, orgs::Array{Organisms.Organism, 1}, t::Int, settings::Dict{String, Any}, orgsref, landpars::Any, tdist::Any)#Setworld.LandPars)}

    lost = Int64[]

    # Only seeds that have been released can disperse 
    for d in seedsi
        if orgs[d].kernel == "short"
            µ, λ = [µ_short λ_short]
        elseif orgs[d].kernel == "medium"
            µ, λ = [µ_medium λ_medium]           
        elseif orgs[d].kernel == "long"
            µ, λ = [µ_long λ_long]           
        elseif orgs[d].kernel in ["medium-short", "short-medium"]
            µ, λ = rand([[µ_short λ_short],
		   	 [µ_medium λ_medium]])
        elseif orgs[d].kernel in ["medium-long", "long-medium"]
            µ, λ = rand([[µ_long λ_long],
		   	 [µ_medium λ_medium]])    
        elseif orgs[d].kernel in ["long-short", "short-long"]
            µ, λ = rand([[µ_short λ_short],
		   	 [µ_long λ_long]])                
        elseif orgs[d].kernel == "any"
            µ,λ = rand([[µ_short λ_short],
	    	        [µ_medium λ_medium],
	    		[µ_long λ_long]])
        else
            error("Check dispersal kernel input for species $(orgs[d].sp).")
        end

	dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ,λ),1)[1])

        # Find the cell to which it is dispersing
        θ = rand(Distributions.Uniform(0,2),1)[1]*pi
        xdest = orgs[d].location[1] + dist*round(Int64, cos(θ), RoundNearestTiesAway)
        ydest = orgs[d].location[2] + dist*round(Int64, sin(θ), RoundNearestTiesAway) 

        # Check if individual can have a chance of establishing there
        if checkbounds(Bool, landavail, xdest, ydest) && landavail[xdest, ydest] == true # checking the suitability first would make more sense but cant be done if cell is out of bounds 

            orgs[d].location = (xdest,ydest) 

        else # if the new location is in an unavailable habitat or outside the landscape, the seed dies

            push!(lost,d)

        end
    end
    #unity test
    #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
    println("Lost $(length(lost))")
    #end
    deleteat!(orgs,lost)
end

"""
                    germinate(org)
                    Seeds have a probability of germinating (`gprob`).
                    """

function germinate(org::Organisms.Organism, T::Float64, settings::Dict{String, Any})
    Bg = org.b0germ * (org.mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
    gprob = 1 - exp(-Bg)

    if gprob < 0
        error("gprob < 0")
    elseif gprob > 1
        error("gprob > 1")
    end

    germ = false
    if 1 == rand(Distributions.Bernoulli(gprob))
        germ = true
    end
    return germ
end

"""
                                    establish!
                    Seed that have already been released (in the current time step, or previously - this is why `seedsi` does not limit who get to establish) and did not die during dispersal can establish.# only after release seed can establish. Part of the establishment actually accounts for the seed falling in an available cell. This is done in the dispersal() function, to avoid computing this function for individuals that should die anyway. When they land in such place, they have a chance of germinating (become seedlings - `j` - simulated by `germinate!`). Seeds that don't germinate stay in the seedbank, while the ones that are older than one year are eliminated.
                                    """
function establish!(orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any}, orgsref, T::Float64)
    #REFERENCE: May et al. 2009
    establishing = find(x -> x.stage == "e" && rem(t,52) > x.seedon, orgs)

    lost = Int64[]
    germinated = 0

    for o in establishing

        if germinate(orgs[o],T,settings)
            orgs[o].stage = "j"
	    open(abspath(joinpath(settings["outputat"],settings["simID"],"eventslog.txt")),"a") do sim
            writedlm(sim, hcat(t, "germination", orgs[o].stage, orgs[o].age))
        end

        end
    end
    
end

"""
                                    survive!(orgs, nogrowth,landscape)
                                    Organism survival depends on total biomass, according to MTE rate. However, the proportionality constants (b_0) used depend on the cause of mortality: competition-related, where
                                    plants in nogrwth are subjected to two probability rates
                                    """
function survive!(orgs::Array{Organisms.Organism,1}, t::Int, cK::Float64, K::Float64, settings::Dict{String, Any}, orgsref, landavail::BitArray{2},T, nogrowth::Array{Int64,1})
    #open(string("EDoutputs/",settings["simID"],"/simulog.txt"), "a") do sim
    #    println(sim,"Running mortality")
    #end

    deaths = Int64[]
    seeds = find(x -> x.stage == "e", orgs)
    mprob = 0
    Bm = 0

    # Density-independent mortality
    for o in 1:length(orgs)

        if sum(values(orgs[o].mass)) <= 0 || orgs[o] in nogrowth #probably unnecessary to verify negative or null weights
            mprob = 1

        elseif o in seeds
                #if orgs[o].age >= orgs[o].bankduration
                #    mprob = 1
            if (rem(t,52) > 38 || rem(t,52) < 12) && # seeds can only die during winter
	       rem(t,52) > orgs[o].seedoff #seeds that are still in the mother plant cant die. If their release season is over, it is certain thatthey are not anymore, even if they have not germinated 
                Bm = orgs[o].b0mort * (orgs[o].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
                #println("Bm: $Bm, b0mort = $(orgs[o].b0mort), seed mass = $(orgs[o].mass["veg"])")
                mprob = 1 - exp(-Bm)
            else
                mprob = 0
            end

    elseif (orgs[o].stage == "a" && orgs[o].age >= orgs[o].span) #oldies die
        mprob = 1
        
    else #calculate mortality for juveniles or adults
        Bm = orgs[o].b0mort *  (orgs[o].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
        mprob = 1 - exp(-Bm)
    end

    # Check mortality rate to probability conversion
    if mprob < 0
        error("mprob < 0")
        #mprob = 0
    elseif mprob > 1            #mprob = 1
        error("mprob > 1")
    end

    if 1 == rand(Distributions.Bernoulli(mprob))
        push!(deaths, o)
        #println("$(orgs[o].stage) dying INDEP.")
	# test
	open(abspath(joinpath(settings["outputat"],settings["simID"],"eventslog.txt")),"a") do sim
            writedlm(sim, hcat(t, "death", orgs[o].stage, orgs[o].age))
        end

    else
        orgs[o].age += 1
    end
end

deleteat!(orgs, deaths) #delete the ones that are already dying due to mortality rate, so that they can´t also die due to density-dependent

open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"a") do sim
    println(sim, "Density-independent mortality: ", length(deaths))
end

deaths = Int64[] # reset before calculating density-dependent mortality

## Density-dependent mortality

if sum(vcat(map(x -> x.mass["veg"], orgs), 0.00001)) > K
    
    locs = map(x -> x.location,orgs) #probably optimizable
    l = DataFrame(locs)
    
    # separate location coordinates and find all individuals that are in the same location as others (by compaing their locations with nonunique(only possible row-wise, not between tuples. This is the only way to get their indexes
    fullcells = find(nonunique(l)) # indexes in l are the same as in orgs, right?
    #open(string("EDoutputs/",settings["simID"],"/simulog.txt"), "a") do sim
    #    println(sim,"# of cell with competition: $(length(fullcells))")
    #end
    if length(fullcells) > 0
        for c in fullcells

            #find plant that are in the same grid
            samegrid = filter(x -> x.location == locs[c], orgs)

            # unity test
            #open(string("EDoutputs/",settings["simID"],"/simulog.txt"), "a") do sim
            println("N of inds in same cell: $(length(samegrid))")
            #end
            #println("Inds in same cell $(locs[c]) : $samegrid")

            # sum their weight to see if > than carrying capacity.
            println("Same grid mass", sum(map(x -> x.mass["veg"], samegrid)))
            # while the sum of weights in a "shared grid" is higher than cell carrying capacity
            while sum(vcat(map(x -> x.mass["veg"],samegrid),0.00001)) > cK #vcat is necessary for avoid that sum throws an error when empty

                # any individual can die
                d = rand(samegrid,1)[1]

                if d.stage == "e"
                    Bm = d.b0mort * (sum(collect(values(d.mass["veg"]))))^(-1/4)*exp(-aE/(Boltz*T))
                elseif d.stage == "j"
                    Bm = d.b0mort * (sum(collect(values(d.mass["veg"]))))^(-1/4)*exp(-aE/(Boltz*T))
                else
                    Bm = d.b0mort * (sum(collect(values(d.mass["veg"]))))^(-1/4)*exp(-aE/(Boltz*T)) 
                end
                # calculate probability
                mprob = 1 - exp(-Bm)
                
                #unity test: Check mortality rate to probability conversion
                if mprob < 0
                    error("mprob < 0")
                    mprob = 0
                elseif mprob > 1
                    mprob = 1
                    error("mprob > 1")
                end

                o = find(x -> x.id == d.id, orgs)

                if 1 == rand(Distributions.Bernoulli(mprob))
                    #currentk -= sum(values(orgs[o].mass)) #update currentk to calculate density-dependent mortality
                    push!(deaths, o)
                    #println("$(orgs[o].stage) dying INDEP.")

                    # test
                    open(abspath(joinpath(settings["outputat"],settings["simID"],"eventslog.txt")),"a") do sim
                        writedlm(sim, hcat(t, "death-K", orgs[o].stage, orgs[o].age))
                    end
		    
                else
                    orgs[o].age += 1
                end

                # update the list of plants sharing grids, to avoid re-caluting mortality for them 
                samegrid = filter(x -> x.location == locs[c], orgs)

            end
        end
    end
    
    deleteat!(orgs, deaths)
end

open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")), "a") do sim
    println(sim,"$(length(deaths)) dying (density-dependent).","\n",
            "Total individuals: ", length(orgs))
end

end

"""
                                    shedd!()
                                    Plants loose their reproductive biomasses at the end of the reproductive season.
                                    """
function shedd!(orgs::Array{Organisms.Organism,1}, orgsref, t::Int)

    flowering = find(x -> (x.mass["repr"] > 0), orgs) #indexing a string returns a Char type, not String. Therefore, p must be Char ('').
    
    for f in flowering
        if rem(t,52) > orgs[f].floroff
            orgs[f].mass["repr"] = 0
        end
    end
end

"""
                    destroyorgs!(orgs)
                    Kill organisms that where in the lost habitat cells.
                    """
function destroyorgs!(orgs::Array{Organisms.Organism,1}, landavail::BitArray{2}, settings::Dict{String,Any})
    
    kills = []

    for o in 1:length(orgs)
        
        if landavail[orgs[o].location...] == false 
            push!(kills,o)
        end        
    end
    
    if length(kills) > 0 # trying to delete at index 0 generates an error
        deleteat!(orgs, kills)
    end
    #unity test
    #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
    #println("Killed orgs: $(length(kills))")
    #end
end

"""
                                    occupied()
                                    Finds indexes of occupied cells in the landscape.
                                    """
function occupied(orgs)
    locs = fill!(Array{Tuple{Int64,Int64,Int64}}(reverse(size(orgs))),(0,0,0))
    #masses = zeros(Float64,size(orgs))
    map!(x -> x.location,locs,orgs) #probably optimizable
    l = DataFrame(locs)
    #println("$l")
    #map!(x -> x.mass["veg"],masses,orgs)
    # separate location coordinates and find all individuals that are in the same location as others (by compaing their locations with nonunique(only possible row-wise, not between tuples. This is the only way to get their indexes
    fullcells = find(nonunique(l)) # indexes in l are the same as in orgs, right?
    #open(string("EDoutputs/",settings["simID"],"/simulog.txt"), "a") do sim
    #    println(sim,"# of cell with competition: $(length(fullcells))")
    #end
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
# function mate!(org::Organism)
#     x, y, frag = org.location #another org of same sp should match the locations of focus
#     sp = org.sp
#
#     # 1. check in the location field of orgs array:
#     # 1.a inside same frag, look for locaions inside the squared area.
#     # 2. when matching, differentiate between autotrphsa and the rest
#     if org.stage == "a"
#
#         for o in 1:length(orgs) #look for partners
#             #TODO optimize indexation of field location in arrray
#             #TODO memory-wise, is it better to put all ifs together?
#             # 1:1 sex-ratio,
#             if frag == orgs[o].location[3]
#                 if orgs[o].location[1:2] in collect(Iterators.product(x-1:x+1,y-1:y+1))
#                     #check sp, self and already reproduced
#                     # if (sp == orgs[o].sp && !(Base.isequal(org, orgs[o])) && org.repr == false && orgs[o].repr = false)
#                     #     #TODO add stochasticity
#                     #     org.repr = true
#                     #     orgs[o].repr = true
#                     #
#                     #     parents_genes = [org.genotype, orgs[o].repr]
#                     # end
#                 end
#             end
#         end
#
#     else
#         continue
#     end
#     return parents_genes #TODO check if it conflicts with modifying orgs
# end


end
