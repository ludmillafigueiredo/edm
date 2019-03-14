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

export Organism, OrgsRef, initorgs, develop!, allocate!, mate!, mkoffspring!, disperse!, germinate, establish!, survive!, shedd!, destroyorgs!, release!

# Set up model constants
const Boltz = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
const aE = 0.63 #0.65 - eV Brown & Sibly MTE book chap 2
const µ_wind = 0.1
const λ_wind = 3
const µ_ant = 1
const λ_ant = 0.2

# Initial organisms parametrization is read from an input file and stored in OrgsRef
mutable struct OrgsRef
    sp_id::Array{String, 1}
    abund::Dict{String,Int}
    kernel::Dict{String,String}
    clonal::Dict{String,Bool}
    emass::Dict{String,Float64}
    elong_mean::Dict{String,Int}
    elong_sd::Dict{String,Float64}
    floron_mean::Dict{String,Int}
    floron_sd::Dict{String,Float64}
    floroff_mean::Dict{String,Int}
    floroff_sd::Dict{String,Float64}
    seedon_mean::Dict{String,Int}
    seedon_sd::Dict{String,Float64}
    seedoff_mean::Dict{String,Int}
    seedoff_sd::Dict{String,Float64}
    span_mean::Dict{String,Int}
    span_sd::Dict{String,Float64}
    b0grow_mean::Dict{String,Float64}
    b0grow_sd::Dict{String,Float64}
    b0m_mean::Dict{String,Float64}
    b0m_sd::Dict{String,Float64}
    b0germ_mean::Dict{String,Float64}
    b0germ_sd::Dict{String,Float64}
    maxseedn_mean::Dict{String,Int}
    maxseedn_sd::Dict{String,Float64}
    maxmass::Dict{String,Int}
end

mutable struct Organism
    id::String
    stage::String #e,j,a
    location::Tuple # (x,y)
    sp::String #sp id, easier to read
    kernel::String
    clonal::Bool
    #### Evolvable traits ####
    emass::Float64
    elong::Int64
    floron::Int
    floroff::Int
    seedon::Int
    seedoff::Int
    span::Int64
    b0grow::Float64
    b0m::Float64
    b0germ::Float64
    wseedn::Int64 # maximum number of seeds produced/week
    maxmass::Float64
    #### State variables #### 
    age::Int64 # control death when older than max. lifespan
    mass::Dict
    mated::Bool
end

"""
                    initorgs(landavail, orgsref,id_counter)

                    Initializes the organisms characterized in the input info stored in `orgsref` and distributes them in the available landscape `landavail`. Stores theindividuals in the `orgs` array, which holds all organisms being simulated at any given time.

                    """
function initorgs(landavail::BitArray{N} where N, orgsref::Organisms.OrgsRef, id_counter::Int)

    orgs = Organism[]
    
    for s in orgsref.sp_id # all fragments are populated from the same species pool

        # create random locations
        XYs = hcat(rand(1:size(landavail,1),orgsref.abund[s]),
                   rand(1:size(landavail,2),orgsref.abund[s]))

        for i in 1:orgsref.abund[s]

            id_counter += 1 # update individual counter
	    non0sd = 1e-7 #avoid error due to sd = 0

            neworg = Organism(hex(id_counter),
                              rand(["a" "j" "e"]),
                              (XYs[i,1],XYs[i,2]),
                              s,
                              orgsref.kernel[s],
                              orgsref.clonal[s],
                              orgsref.emass[s],
                              Int(round(rand(Distributions.Normal(orgsref.elong_mean[s],orgsref.elong_sd[s]+non0sd),1)[1])),
                              Int(round(rand(Distributions.Normal(orgsref.floron_mean[s],orgsref.floron_sd[s]+non0sd),1)[1])),
                              Int(round(rand(Distributions.Normal(orgsref.floroff_mean[s],orgsref.floroff_sd[s]+non0sd),1)[1])),
                              Int(round(rand(Distributions.Normal(orgsref.seedon_mean[s],orgsref.seedon_sd[s]+non0sd),1)[1])),
                              Int(round(rand(Distributions.Normal(orgsref.seedoff_mean[s],orgsref.seedoff_sd[s]+non0sd),1)[1])),
                              Int(round(rand(Distributions.Normal(orgsref.span_mean[s], orgsref.span_sd[s]+non0sd),1)[1])),
                              Int(round(rand(Distributions.Normal(orgsref.b0grow_mean[s],orgsref.b0grow_sd[s]+non0sd),1)[1])),
                              Int(round(rand(Distributions.Normal(orgsref.b0m_mean[s],orgsref.b0m_sd[s]+non0sd),1)[1])),
                              Int(round(rand(Distributions.Normal(orgsref.b0germ_mean[s],orgsref.b0germ_sd[s]+non0sd),1)[1])),
                              0, #wseedn
                              0.0, #maxmass
                              0, #age
                              Dict("veg" => 0.0, "repr" => 0.0), #mass
                              false) #mated

            ## Set conditional traits and variables
            # weekly number of seeds
	    maxseedn = Int(round(rand(Distributions.Normal(orgsref.seedoff_mean[s],orgsref.seedoff_sd[s]+non0sd),1)[1]))
            neworg.wseedn = Int(round(maxseedn/(neworg.floroff-neworg.floron+1)))

            # initial biomass
            if neworg.stage == "e"                  
                neworg.mass["veg"] = neworg.emass
                neworg.age = 1
            elseif neworg.stage in ["j" "a"]
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
function allocate!(orgs::Array{Organism,1}, t::Int64, aE::Float64, Boltz::Float64, settings::Dict{String, Any},orgsref::Organisms.OrgsRef,T::Float64)
    #1. Initialize storage of those that dont growi and will have higher prob of dying (later)
    nogrowth = Int64[]

    growing = find(x->(x.stage == "a" || x.stage == "j"),orgs)
    for o in growing

        b0 = orgs[o].b0grow
        
        #only vegetative biomass helps growth
        grown_mass = b0*(orgs[o].mass["veg"])^(3/4)*exp(-aE/(Boltz*T))

        if isapprox(grown_mass,0) # if it is not growing, there is no need to allocate
            push!(nogrowth,o)
        elseif orgs[o].stage == "j"
            # juveniles grow vegetative biomass only
            orgs[o].mass["veg"] += grown_mass 
        elseif orgs[o].stage == "a" &&
            (orgs[o].floron <= rem(t,52) < orgs[o].floroff) &&
            (sum(collect(values(orgs[o].mass))) >= 0.5*(orgs[o].maxmass))
            # adults in their reproductive season and with enough weight, invest in reproduction
            #sowingmass = (5.5*(10.0^(-2)))*((orgs[o].mass["veg"]/(orgs[o].floroff-orgs[o].floron + 1) + grown_mass)^0.95)
            if haskey(orgs[o].mass,"repr")
                orgs[o].mass["repr"] += grown_mass #sowingmass 
            else
                orgs[o].mass["repr"] = grown_mass #sowingmass
            end
        elseif orgs[o].stage == "a" && orgs[o].mass["veg"] < orgs[o].maxmass
            # adults that have not yet reached maximum size can still grow vegetative biomass, independently of the season
            orgs[o].mass["veg"] += grown_mass 
        end            
    end
    return nogrowth
end

"""
                develop!()
                Controls individual juvenile maturation.
                """
function develop!(orgs::Array{Organism,1}, orgsref::Organisms.OrgsRef)
    juvs = find(x->x.stage == "j",orgs)

    for j in juvs
        if  orgs[j].mass["veg"] >= 0.5*orgs[j].maxmass
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
function mate!(orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any}, scen::String, regime:: String, td::Array{Int64,1}, remaining)

    ready = find(x-> x.stage == "a" && x.mass["repr"] > x.emass, orgs) # TODO find those with higher reproductive mas than the mean nb of seeds * seed mass.
    pollinated = []
    npoll = 0

    if length(ready) > 0 # check if there is anyone flowering
        # Scenarios were pollination is not species-specific
        if scen in ["indep" "equal"]  # calculation of number of pollinated individuals is different, but the actual pollination (non species-specific) is

            # Base number of pollinated flowes
            defaultnpoll = rand(Distributions.Binomial(Int(ceil(length(ready) * 0.01)),0.5))[1] 
            
            # Determine number of individuals that get pollinated (species is not relevant)
            if scen == "indep"
                npoll = defaultnpoll # Fishman & Hadany's proportion of visited flowers
                println("Scenario of INDEP pollination loss")

            elseif scen == "equal" #all species lose pollination randomly (not species-specific)
                println("Scenario of EQUAL pollination loss")
                # Check the regime
                if t in td 
                    if regime == "pulse"
                        npoll = Int(ceil(defaultnpoll * remaining))
                    elseif regime == "exp"
                        npoll = Int(ceil(defaultnpoll * exp(-0.5*(t-td))))
                        # exp(tp - t) makes the pollination loss decrease from 1 (tp = t) to 0
                        println("Regime of EXP pollination loss, with $(length(ready)) being ready and $npoll being pollinated")
                    elseif regime == "press"
                        npoll = Int(ceil(defaultnpoll * remaining[find(td.==t)]))
                    else
                        error("Please choose a pollination regime \"regime\" in insect.csv:
                                - \"pulse\":
                                - \"press\":
                                - \"exp\"")                    
                    end
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
function mkoffspring!(orgs::Array{Organisms.Organism,1}, t::Int64, settings::Dict{String, Any},orgsref::Organisms.OrgsRef, id_counter::Int, landavail::BitArray{2})

    offspring = Organism[]

    # Asexually produced offspring
    asexuals = filter(x -> x.mated == false && x.clonal == true && x.mass["repr"] > x.emass, orgs) #find which the individuals the can reproduce assexually and then go through them, by species

    for sp in unique(getfield.(asexuals, :sp))
        cloning = find(x -> x.sp == sp && x.id in unique(getfield.(asexuals, :id)) , orgs)

        # start production counting
        spclonescounter = 0

        for c in cloning # mothers cloning
            offs = div(0.5*orgs[c].mass["repr"], orgs[c].emass)

            # unity test
            if  orgs[c].mass["repr"] <= 0
                error("Negative reproductive biomass") #because offs is an integer, reproductive biomass should not become negative
            end

            if offs <= 0
                continue
            else
                # limit offspring production to the maximal number of seeds the species can produce
                offs > orgs[c].wseedn ? offs = orgs[c].wseedn : offs

                # count clonal production  
                spclonescounter += offs

                # update reproductive mass
                orgs[c].mass["repr"] -= (offs * orgs[c].emass)

                # get a copy of the mother, which the clones will look like 
                clonetemplate = deepcopy(orgs[c])
                clonetemplate.stage = "j" #clones have already germinated
                clonetemplate.mass["veg"] = orgs[c].emass
                clonetemplate.mass["repr"] = 0.0

                for o in offs

                    clone = deepcopy(clonetemplate)
                    # check if the new location is actually available before creating the clone
                    clone.location = (clonetemplate.location[1] + rand(Distributions.Bernoulli())[1],
                                      clonetemplate.location[2] + rand(Distributions.Bernoulli())[1]) # clones are spread in one of the neighboring cells - or in the same as the mother

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

    # Sexually produced offspring
    ferts = filter(x -> x.mated == true, orgs)

    for sp in unique(getfield.(ferts, :sp))

        sowing = find(x -> x.sp == sp || x.id in unique(getfield.(ferts, :id)), orgs)

        spoffspringcounter = 0

        for s in sowing

            emu = orgs[s].emass
            offs = div(0.5*orgs[s].mass["repr"], emu)

            if offs <= 0
                continue
            else
                # limit offspring production to the maximal number of seeds the species can produce
                offs > orgs[s].wseedn ? offs = orgs[s].wseedn : offs

                # count species offspring for output
                spoffspringcounter += offs

                orgs[s].mass["repr"] -= (offs * emu)

                # unity test
                if  orgs[s].mass["repr"] <= 0
                    error("Negative reproductive biomass")
                end

                # get another parent
                conspp = orgs[rand(find(x -> x.sp == sp && x.stage == "a", orgs))] 

                for n in 1:offs

                    id_counter += 1 # update individual counter

                    embryo = deepcopy(orgs[s])

                    #newvalue = rand(Distributions.Normal(0,abs(embryo.emass-conspp.emass)/embryo.emass))
                    #embryo.emass + newvalue >= orgsref.emass[embryo.sp] ? # if seed biomass or minimal biomass would smaller than zero, it does not chenge
                    #embryo.emass += newvalue : embryo.emass += 0
                    embryo.b0grow += rand(Distributions.Normal(0,abs(embryo.b0grow-conspp.b0grow+0.0000001)/embryo.b0grow))
                    embryo.b0m += rand(Distributions.Normal(0,abs(embryo.b0m-conspp.b0m+0.0000001)/embryo.b0m))
                    embryo.b0germ += rand(Distributions.Normal(0,abs(embryo.b0germ-conspp.b0germ+0.0000001)/embryo.b0germ))
                    embryo.floron += Int(round(rand(Distributions.Normal(0,abs(embryo.floron-conspp.floron+0.0000001)/embryo.floron)),RoundUp))
                    embryo.floroff += Int(round(rand(Distributions.Normal(0,abs(embryo.floroff-conspp.floroff+0.0000001)/conspp.floroff)),RoundUp))
                    embryo.seedon += Int(round(rand(Distributions.Normal(0,abs(embryo.seedon-conspp.seedon+0.0000001)/embryo.seedon)),RoundUp))
                    embryo.seedoff += Int(round(rand(Distributions.Normal(0,abs(embryo.seedoff-conspp.seedoff+0.0000001)/embryo.seedoff)),RoundUp))
                    
                    # set embryos own individual non-evolvable traits
                    embryo.id = hex(id_counter) 
                    embryo.mass = Dict("veg" => embryo.emass,
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
            writedlm(seedfile, hcat(t, sp, "s", "sex", spoffspringcounter))
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
function release!(orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any},orgsref::Organisms.OrgsRef )

    # Individuals being released in any given week are: in embryo stage (=seed= & in their seed release period (seedon <= t <= seedoff for the species)

    seedsi = find(x -> x.stage == "e" && x.age == 0 && x.seedon <= rem(t,52) < x.seedoff, orgs)
    # using a condition "outside" orgs might not work. This condition with orgsref only works because orgsref always has the sp names of x as keys in the dictionnary. If presented with a key that it does do contain, it throws an error.

    return seedsi
end

"""
                disperse!(landscape, orgs, t, seetings, orgsref,)
                Seeds are dispersed.
                """

function disperse!(landavail::BitArray{2}, seedsi, orgs::Array{Organisms.Organism, 1}, t::Int, settings::Dict{String, Any}, orgsref::Organisms.OrgsRef, landpars::Any, tdist::Any)#Setworld.LandPars)}

    lost = Int64[]

    # Only seeds that have been released can disperse 
    for d in seedsi
        if orgs[d].kernel == "a"
            dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ_ant,λ_ant),1)[1])
        elseif orgs[d].kernel == "w"
            dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ_wind,λ_wind),1)[1])           
        elseif orgs[d].kernel == "wa"
            µ,λ = rand([[µ_wind λ_wind],[µ_ant λ_ant]])
            dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ,λ),1)[1])

        else
            error("Check dispersal kernel input for species $(orgs[d].sp).")
        end

        # Find the cell to which it is dispersing
        θ = rand([0 0.5π π 1.5π]) # TODO tentar rand(0:2)*pi?
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
function establish!(orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any}, orgsref::Organisms.OrgsRef, T::Float64)
    #REFERENCE: May et al. 2009
    establishing = find(x -> x.stage == "e" && rem(t,52) > x.seedon, orgs)

    lost = Int64[]

    for o in establishing

        if germinate(orgs[o],T,settings)
            orgs[o].stage = "j"
        end

    end
end

"""
                                survive!(orgs, nogrowth,landscape)
                                Organism survival depends on total biomass, according to MTE rate. However, the proportionality constants (b_0) used depend on the cause of mortality: competition-related, where
                                plants in nogrwth are subjected to two probability rates
                                """
function survive!(orgs::Array{Organisms.Organism,1}, t::Int, cK::Float64, K::Float64, settings::Dict{String, Any}, orgsref::Organisms.OrgsRef, landavail::BitArray{2},T, nogrowth::Array{Int64,1})
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
            if orgs[o].age >= orgs[o].elong
                mprob = 1
            elseif rem(t,52) > orgs[o].seedoff #seeds that are still in the mother plant cant die. If their release season is over, it is certain thatthey are not anymore, even if they have not germinated 
                Bm = orgs[o].b0m * (orgs[o].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
                #println("Bm: $Bm, b0m = $(orgs[o].b0m), seed mass = $(orgs[o].mass["veg"])")
                mprob = 1 - exp(-Bm)
            else
                mprob = 0
            end 

        elseif orgs[o].age >= orgs[o].span #oldies die
            mprob = 1

        else
            Bm = orgs[o].b0m * (orgs[o].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
            mprob = 1 - exp(-Bm)                         
        end

        # Check mortality rate to probability conversion
        if mprob < 0
            error("mprob < 0")
            #mprob = 0
        elseif mprob > 1
            #mprob = 1
            error("mprob > 1")
        end

        if 1 == rand(Distributions.Bernoulli(mprob))
            push!(deaths, o)
            #println("$(orgs[o].stage) dying INDEP.")
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

                    if d.stage == "a"
                        Bm = d.b0m * (sum(collect(values(d.mass["veg"]))))^(-1/4)*exp(-aE/(Boltz*T))
                    elseif d.stage == "j"
                        Bm = d.b0m * (sum(collect(values(d.mass["veg"]))))^(-1/4)*exp(-aE/(Boltz*T))
                    else
                        Bm = d.b0m * (sum(collect(values(d.mass["veg"]))))^(-1/4)*exp(-aE/(Boltz*T))
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
function shedd!(orgs::Array{Organisms.Organism,1}, orgsref::Organisms.OrgsRef, t::Int)

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
