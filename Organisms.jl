"""
This module contains the

Organisms have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the range of search).
"""
module Organisms

using Distributions
using Setworld
using Fileprep

export Organism, InitOrgs, newOrgs, projvegmass!, compete, develop!, allocate!, meanExP, checkboundaries, reproduce!, disperse!, germinate, establish!, survive!

#TODO put them in OrgsRef
const Boltz = 8.62e-5 # eV/K Brown & Sibly MTE book chap 2
const aE = 0.65 # eV Brown & Sibly MTE book chap 2
const plants_gb0 = (10^(10.15))/40 # 10e10.15 is the annual plant biomass production (Ernest et al. 2003) transformed to weekly base, with growth not happening during winter
const plants_mb0 = 1.5029220413821088e11 #adjustted accordung to 1 death per individual for 1g (MTEpar notebook)
const plants_fb0 = exp(30.0) # fertility rate
const seedmassµ = 0.8

# Initial organisms parametrization
mutable struct InitOrgs
    fgroups::Any
    sps::Any
    init_stage::Any
    init_abund::Any
    #genotypes #TODO initialie those in functions
    biomassμ::Any
    biomasssd::Any
    dispμ::Any
    dispshp::Any
    radius::Any
    InitOrgs() = new() #TODO check if new() is necessary
end

mutable struct Organism
    id::String
    location::Tuple # (x,y,frag)
    sp::String
    stage::String #e,j,a
    age::Int64 # controls passing stages, phenology and
    reped::Bool # reproduced?
    fgroup::String # Plant, insect or more precise functional group
    genotype::Array{String,2}
    biomass::Dict
    disp::Tuple #dispersal kernel parameters (a and b) TODO tuple
    radius::Int64 # TODO reproductive and vegetative area of influence. Not Tuple because not
    #Organism() = new()
end

"""
    newOrg(fgroups, init_abund, biomassμ, biomasssd)
newOrg() creates new `init_abund` individuals of each  functional group (`fgroups`), with biomass (N(`biomassμ`,`biomasssd`)). It is called during the initialization and all  the Organisms parameters are read from a parameter file, for the specific functional group. All organisms, from all funcitonal groups, are stored in the same array `orgs` thathte function returns.

    newOrg(parents, reprmode)
 newOrg() generates offspringn, after reproduc- newOrg() generates offspringn, after reproduction. For offspring creation, the `genotype` field is filled according to parental genotype and mode of reproduction and the other fields are copied from the parents.
# Methods
1. Upon initialization: `newOrg(fgroups, init_abund, biomassμ, landscape)`
`fgroups` list of funcitonal groups to be simulated
`init_abund` initial abundances
`biomassμ` mean biomass

# Arguments:
`reprmode::String` specificies how the new organisms are generated:
    `init` initialization of the model
    `sex` sexual reproduction
    `clone` assexual reproduction
`parent_s::Array{Organism,N}` array with single parent for clones, both for sexual reproduction
`quant::Int64` is nb of new individuals or offspring to be created
"""
function newOrgs(landscape::Array{Setworld.WorldCell,3},initorgs::Organisms.InitOrgs)

    orgs = Organism[]

    for frag in 1:size(landscape,3)

        for f in 1:length(initorgs.fgroups) #TODO check the organisms file format: so far, all fragments get the same sps

            XYs = hcat(rand(1:size(landscape,1),initorgs.init_abund[f]),
                       rand(1:size(landscape,2),initorgs.init_abund[f]))

            for i in 1:initorgs.init_abund[f]
                neworg = Organism(string(initorgs.fgroups[f][1:3], length(orgs) + 1),
                                  (XYs[i,1],XYs[i,2],frag),
                                  initorgs.sps[f], #TODO in order to have redundant species, this should take a rdm? sps inside the sps pool of that fragment OR o input the sps teria que ser casado com o de funcitonal group e abundance (cada linha eh uma sp, com abundancia etc...)
                                  initorgs.init_stage[f],
                                  0,
                                  false,
                                  initorgs.fgroups[f],
                                  ["" ""], #initialize with function
                                  Dict("veg" => rand(Distributions.Normal(initorgs.biomassμ[f],initorgs.biomasssd[f]))),
                                  (initorgs.dispμ[f], initorgs.dispshp[f]),
                                  initorgs.radius[f])

                push!(orgs, neworg)

                #global IDcounter += 1
            end
        end
    end
    #TODO different fgroups have different metabolic rates. Get it form OrgsRef
    #TODO or make sure that it has varied in the global scope as well rather export it?
   return orgs
end

# function newOrgs(reprmode::String, parent_s::Array{Organism, N} where N, quant:Int64)
# TODO integrate it with reproduce!()
#     #     # TODO arguments are quite different. Should it be a different function?
#     #     # # TODO For the reproduction method of this fct: The genotype and individual location depends on the group of reproduction
#     for i in 1:quant
#         neworg = Organism(string(parent_s[1].fgroup[1:3], length(orgs) + 1),
#                           (XYs[i,1] XYs[i,2] frag),
#                           initorgs.sps[f],
#                           initorgs.init_stage[f],
#                           0,
#                           false,
#                           initorgs.fgroups[f],
#                           ["placeholder" "placeholder"], #initialize with function
#                           Dict("veg" => rand(Distributions.Normal(initorgs.biomassμ[f],initorgs.biomasssd[f]))),
#                           (initorgs.dispμ[f] initorgs.dispshp[f]),
#                           [initorgs.radius[f]])
#         push!(orgs, neworg)
#
#         #global IDcounter += 1
#     end
# end

"""
    projvegmass!(landscape,orgs)
Rewrites the projected mass of each organisms stored in `orgs` into the `neighs` field of `landscape`. This projection means that the total biomass is divided into the square area delimited by the organism's `radius`.
"""
function projvegmass!(landscape::Array{Setworld.WorldCell, 3}, orgs::Array{Organism,N} where N)
    # empt neighs to rewrite
    # TODO more efficient?
    for f in 1:size(landscape,3)
        for y in 1:size(landscape,2), x in 1:size(landscape,1)
            if length(keys(landscape[x,y,f].neighs)) == 0
                continue
            else
                landscape[x,y,f].neighs = Dict()
            end
        end
    end

    for o in 1:length(orgs)
        x, y, frag = orgs[o].location
        orgs[o].radius = round(Int64, (sqrt(orgs[o].biomass["veg"]^(2/3)) - 1)/2, RoundNearestTiesAway)
        r = orgs[o].radius # separated for debugging
        # unity test
        println("radius $r")
        fg = orgs[o].fgroup
        projmass = /(orgs[o].biomass["veg"], ((2*r+1)^2))

        #unity test
        #println(simulog, orgs[o].id, " has biomass", orgs[o].biomass, " and projects ", projmass)

        for j in (y-r):(y+r), i in (x-r):(x+r) #TODO usar a funcao da FON Q trabalha com quadrantes? dar mais peso para steming point?
            if !checkbounds(Bool,landscape[:,:,frag],j,i) # check boundaries: absorbing borders: the biomass is not re-divided to the amount of cells inside the fragment. What is projected outside the fragmetn is actually lost: Edge effect
                continue
            else
                if haskey(landscape[i,j,frag].neighs,fg)
                    landscape[i,j,frag].neighs[fg] += projmass
                else
                    landscape[i,j,frag].neighs[fg] = projmass
                end
            end
        end
    end
end

"""
    compete(landscape, orgs)
Each plant in `orgs` will check neighboring cells (inside it's zone of influence radius `r`) and detected overlaying ZOIs. When positive, the proportion of 'free' plant biomass is calculated: (focus plant biomass - sum(non-focus vegetative biomass))/(focus plant biomass), and normalized to the total area of projected biomass .
"""
function compete(landscape::Array{Setworld.WorldCell, 3}, org::Organism)

    x, y, frag = org.location
    fg = org.fgroup
    r = org.radius

    compterm = 0
    nbsum = 0

    # 2. Look for neighbors in the square area (2r+1) delimited by r
    for j in (y-r):(y+r), i in (x-r):(x+r)
        if !checkbounds(Bool,landscape[:,:,frag],i,j)
            continue
            #unity test
            #println(simulog, orgs.id, " is projecting outside the landscape, because it is at: ", org[o].location)
        elseif haskey(landscape[i,j,frag].neighs,fg)
            #landscape[i,j,frag].neighs[fg] > 0 # check the neighborhood of same fgroup for competition
            nbsum += /(landscape[i,j,frag].neighs[fg] - org.biomass["veg"],(2*r+1)^2) #sum vegetative biomass of neighbors only (exclude focus plant own biomass)
        else #!haskey(landscape[i,j,frag].neighs,fg) #no competition
            nbsum += 0
        end
    end

    compterm = /(org.biomass["veg"] - nbsum, org.biomass["veg"])
    # unity test
    println("compterm $compterm")

    return compterm
end

"""
    allocate!(orgs, landscape, aE, Boltz, OrgsRef)
Calculates biomass gain according to MTE rate and depending on competition. Competition is measured via a biomass-based index `compterm` (`compete` function). This term gives the proportion of actual biomass gain an individual has. If competition is too strong (`compterm` < 0), the individual has a higher probability of dying. If not, the biomass is allocated to growth or reproduction, according to the developmental `stage` of the organism and the season (`t`) (plants start allocating to week 12).
"""
function allocate!(landscape::Array{Setworld.WorldCell,3}, orgs::Array{Organism,N} where N, t::Int64, aE::Float64, Boltz::Float64)
    #1. Initialize storage of those that are ont growing and have higher prob of dying (later)
    nogrowth = Int64[]

    #2. Calculate growth for all plants (#TODO filter insects out)
    for o in 1:length(orgs)
        # 2.a Get local temperature
        T = landscape[orgs[o].location[1], orgs[o].location[2], orgs[o].location[3]].temp

        # This MTE rate comes from dry weights: fat storage and whatever reproductive structures too, but not maintenance explicitly
        # Any cost related to insufficient minimal biomass goes into the survival probability function
        # 2.b Check for competition
        compterm = compete(landscape, orgs[o])
        # unity test
        #println(simulog, org.id," weights",org.biomass["veg"]," had $nbsum g overlap")
        open("EDoutputs/simulog.txt","a") do sim
            println(sim, orgs[o].id, "  compterm $compterm")
        end

        if compterm < 0
            #2.c
            push!(nogrowth,o)
        else
            #unity test
            println("current biomass ", sum(values(orgs[o].biomass)))

            grown_mass = plants_gb0 * (compterm *sum(values(orgs[o].biomass)))^(3/4) * exp(-aE/(Boltz*T))

            # unity test
            println("gained ",compterm *sum(values(orgs[o].biomass)))

            #Resource allocation schedule
            #TODO make it more ellaborate and includde trade-offs
            if orgs[o].stage == "e"  #TODO accountant for resistant stages?
                # embryos only consume reserves: TODO more realistic
                orgs[o].biomass["veg"] -= 0.01*org[o].biomass["veg"]

            elseif orgs[o].stage == "j"
                # juveniles grow
                orgs[o].biomass["veg"] += grown_mass

            elseif orgs[o].stage == "a" # && 12 <= rem(t, 52) < 25
                # adults reproduc
                if haskey(orgs[o].biomass,"reprd")
                    orgs[o].biomass["reprd"] += grown_mass
                else
                    orgs[o].biomass["reprd"] = grown_mass
                end
            else
                continue
            end
        end
    end
    # unity test
    open("EDoutputs/simulog.txt","a") do sim
        println(sim, "Not growing $nogrowth")
    end
    return nogrowth
end

"""
    develop!()
Controls individual stage transition.
"""
function develop!(orgs::Array{Organism,N} where N)
    #TODO plant growth
    juvs = find(x->x.stage == "j",orgs)

    for j in juvs
        orgs[j].stage = "a"
    end
    #TODO insect holometabolic growth
    #TODO insect hemimetabolic development
end
"""
    meanExP(a,b)
Draws a mean dispersed distance from the Exponential Power dispersal kernel.
"""
function meanExP(a::Float64,b::Float64)
    dist = a*rand(Distributions.Gamma(3/b))/rand(Distributions.Gamma(2/b))
    return dist
end

"""
    checkboundaries(sourcefrag,xdest, ydest, fdest)
`source` and `dest` contain the location indexes of the source (mother plant) and the pollen/seed. `checkboundaires()` verifies whether the new polen/seed location `(x,y)` is inside a habitat fragment (same as the source -`frag`- or another one insed the patch). Return a boolean that controls whether the process (reproduction or emergency/germination) proceeds or not.
"""
function checkboundaries(landscape::Array{Setworld.WorldCell,3}, xdest::Int64, ydest::Int64, fdest::Int64)
    #check inside frag
    if checkbounds(Bool, landscape[:,:,fdest], xdest, ydest)
    #if (dx <= size(landscape[:,:,sf])[1] && dy <= size(landscape[:,:,sf])[2])
        inbound = true
    #elseif (dx <= size(landscape[:,:,sf])[1] && dy <= size(landscape[:,:,sf])[2])
    #TODO check ouside frag (more complex than checkbounds): how to detect the direction of neighboring fragments? use a matrix dependent probability?
        # -> to find fdest: compare dist with matriyes sizes and call possible fragmetns in the third index of location)
        # and checkbounds(landscape[:,:,:], xdest, ydest, fdest)
        #inbound = true
    else
        inbound = false
    end
    return inbound
end

"""
    reproduce!(landscape,orgs)
Assigns proper reproduction mode according to the organism functional group. This controls whether reproduction happens or not, for a given individual: plants depend on pollination, while insects do not. Following, it handles fertilization of new embryos and calculates offspring production. New individuals are included in the community at the end of current timestep.
"""
function reproduce!(landscape::Array{Setworld.WorldCell, 3}, orgs::Array{Organisms.Organism,N} where N, t::Int64)
    #TODO sort out reproduction mode (pollination or not) according to functional group
    # if # pollination depending plants
    #     pollination()
    # elseif "insect"
    #     parents_genes = mate!()
    # end

    reproducing = filter(x -> x.stage == "a" && haskey(x.biomass, "reprd"), orgs)

    #unity test
    open("EDoutputs/simulog.txt","a") do sim
        println(sim, "Reproducing: $reproducing week $t")
    end

    offspring = Organism[]

    for o in 1:length(reproducing)

        T = landscape[reproducing[o].location[1], reproducing[o].location[2], reproducing[o].location[3]].temp

        offsprgB = round(Int64, plants_fb0 * sum(values(reproducing[o].biomass))^(-1/4) * exp(-aE/(Boltz*T)),  RoundNearestTiesAway) #TODO stochasticity

        #unity test
        open("EDoutputs/simulog.txt","a") do sim
            println(sim, "Offspring of ", orgs[o], ": ",offsprgB)
        end

        for n in 1:offsprgB
            #TODO check for a quicker way of creating several objects of composite-type
            embryo = Organism(string(reproducing[o].fgroup[1:3], length(orgs) + length(offspring) + 1),
            reproducing[o].location, #location is given according to functional group and dispersal strategy, in disperse!()
            reproducing[o].sp,
            "e",
            0,
            false,
            reproducing[o].fgroup,
            ["placeholder" "placeholder"], #come from function
            #rand(Distributions.Normal(OrgsRef.seedbiomassμ[reproducing[o].fgroup],OrgsRef.biomasssd[reproducing[o].fgroup])),
            #[OrgsRef.dispμ[f] OrgsRef.dispshp[f]],
            #OrgsRef.radius[f])
            Dict("veg" => reproducing[o].biomass["veg"]*0.01), #TODO use seed size for the fgroup
            reproducing[o].disp,
            reproducing[o].radius) # could be 0, should depend on biomass

            push!(offspring, embryo)
        end

        reproducing[o].biomass["reprd"] -= (offsprgB * seedmassµ)
    end
    append!(orgs, offspring)
    # return offspring
end

"""
    disperse!(offspring)
Butterflies, bees and seeds can/are disperse(d).
"""
function disperse!(landscape::Array{Setworld.WorldCell,3},orgs::Array{Organisms.Organism, N} where N)
    # Seeds (Bullock et al. JEcol 2017)
    # Get seeds from org storage
    dispersing = find(x -> x.stage == "e" && x.age == 0, orgs)
    # takes only seeds not in the seed bank (age>0)
    #TODO include insects
    # Exponential kernel for Ant pollinated herbs: mean = a.(Gamma(3/b)/Gamma(2/b))

    #unity test
    open("EDoutputs/simulog.txt","a") do sim
        println(sim, "Dispersing: $dispersing week $t")
    end

    lost = Int64[]

    for d in dispersing
        # Dispersal distance from kernel:
        if orgs[d].fgroup == "wind"
            #Exp, herbs + appendage #TODO LogSech distribution
            dist = Fileprep.lengthtocell(meanExP(4.7e-5,0.2336))
            #(rand(collect(0.388:0.001:3.226))) Bullock's 50th - 95th percentile
            #higher 9th percentile than ant (when comparing + appendage and 10-36 mg)
            # Check for available habitat (inside or inter fragment)
        elseif orgs[d].fgroup == "ant"
            #ExPherbs, 10-36 mg
            dist = Fileprep.lengthtocell(meanExP(0.3726,1.1615))
            #for tests: (rand(collect(0.499:0.001:1,3056))) Bullock's 50th - 95th percentile
        end

        # Find patch
        θ = rand([0,0.5π,π,1.5π]) #get radian angle of distribution
        xdest = round(Int64, orgs[d].location[1] + dist*sin(θ), RoundNearestTiesAway)
        ydest = round(Int64, orgs[d].location[2] + dist*cos(θ), RoundNearestTiesAway)
        fdest = orgs[d].location[3] #TODO is landing inside the same fragment as the source, for now

        if checkboundaries(landscape, xdest, ydest, fdest)
            orgs[d].location = (xdest,ydest,fdest)
            #unity test
            open("EDoutputs/simulog.txt","a") do sim
                println(sim, org[o].id," dispersed $dist")
            end
        else
            push!(lost,d)
            #unity test
            open("EDoutputs/simulog.txt","a") do sim
                println(sim, org[o].id," dispersed $dist but died")
            end
        end
    end
    deleteat!(orgs,lost)
end

"""
    germinate(org)
Seeds have a probability of germinating (`gprob`).
"""
function germinate(org::Organisms.Organism)
    gprob = 0.5
    germ = false
    if 1 == rand(Distributions.Binomial(1,gprob))
        germ = true
    end
    return germ
end

"""
    establish!
Seeds only have a chance of establishing in patches not already occupied by the same funcitonal group, in. When they land in such place, they have a chance of germinating (become seedlings - `j` - simulated by `germinate!`). Seeds that don't germinate stay in the seedbank, while the ones that are older than one year are eliminated.
"""
function establish!(landscape::Array{Setworld.WorldCell,3}, orgs::Array{Organisms.Organism, N} where N, t::Int64)
    #REFERENCE: May et al. 2009
    #check neighs

    #TODO weekly timing of establishment
    establishing = find(x -> x.stage == "e", orgs)

    #unity test
    open("EDoutputs/simulog.txt","a") do sim
        println(sim, "Establishing seeds: $establishing")
    end

    lost = Int64[]

    for o in establishing
        orgcell = orgs[o].location
        orgfg= orgs[o].fgroup
        if haskey(landscape[orgcell[1], orgcell[2], orgcell[3]].neighs,orgfg)
            push!(lost,o)
        else
            if germinate(orgs[o])
                orgs[o].stage = "j"
                # the ones that dont germinate but are older than 1 year die
            elseif orgs[o].age > 52
                push!(lost,o)
            end
        end
    end
    deleteat!(orgs,lost)
end

"""
    survive!(ors, nogrowth,landscape)
Organism survival depends on total biomass, according to MTE rate. However, the proportionality constants (b_0) used depend on the cause of mortality: competition-related, where
plants in nogrwth are subjected to two probability rates
"""
function survive!(landscape::Array{Setworld.WorldCell,3},orgs::Array{Organisms.Organism,N} where N, nogrowth::Array{Int64,N} where N)

    deaths = Int64[]

    for o in 1:length(orgs)

        if orgs[o].age == 52 # annual plants die
            push!(deaths, o)
        else
            T = landscape[orgs[o].location[1], orgs[o].location[2], orgs[o].location[3]].temp

            mortalconst = plants_mb0 #TODO call it from OrgsRef, when with different functional groups
            mB = mortalconst * (sum(values(orgs[o].biomass)))^(-1/4)*exp(-aE/(Boltz*T))
            mprob = 1 - exp(-mB)

            #unity test
            if rand() > 0.7
                open("EDoutputs/simulog.txt","a") do sim
                    println(sim,orgs[o].id,"-",orgs[o].stage, " mprob: $mprob")
                end
            end

            # individuals that didnt grow have
            if o in 1:length(nogrowth)
                compmortconst = plants_mb0 #TODO use different b_0 for mortality consequence of competition
                cmB = compmortconst * (sum(values(orgs[o].biomass)))^(-1/4)*exp(-aE/Boltz*(T))
                cmprob =  1 - e^(-cmB)
                mprob += cmprob
            end

            if rand(Distributions.Binomial(1,mprob),1) == 1
                #mprob > rand(PoissonBinomial([mprob]),1)[1] # TODO should use a more ellaboratedistribution model? successes are death events, because they are the value that is going to be relavant in here: the amount of individuals to be taken out
                push!(deaths, o)
            else
                orgs[o].age += 1
            end
        end
    end
    deleteat!(orgs, deaths)
    return orgs
end

# """
#     emerge!()
# Germination for seeds, emergency for eggs.
# """
# function emerge!()
#     embryos = filter.( x -> x.stage == "e", orgs)
#     deaths = []
#
#     for o in 1:length(embryos)
#         if embryos[o].fgroup in ["plants"] #TODO sould be more dynamics and read from inputed fgroups
#             germb0 = plants_gb0 #TODO call it from OrgsRe
#         elseif embryos[o].fgroup in ["insects"]
#             germb0 = plants_gb0 #TODO call it from OrgsRe
#         end
#
#         gB = germb0 * (sum(values(embryos[o].biomass)))^(-1/4)*exp(-aE/Boltz*T)
#         gprob = 1 - e^(-gB*embryos[o].age) #TODO this age thing
#
#         if gprob > rand(PoissonBinomial(gprob)) # successes are germinations
#             embryos[o].stage = "j"
#             embryos[o].age = 1
#         else
#             if (embryos[o].fgroup in ["plants"] &&  embryos[o].age <= 52) #TODO seed bank of one year. Make it more complex (count age)
#                 continue
#             else
#                 push!(deaths, o)
#             end
#         end
#     end
#     deleteat!(orgs, deaths)
# end

"""
    pollination!()
Simulates plant-insect encounters and effective pollen transfer.
"""
#TODO find a not too cumbersome way of modelling pollen transfer: store interactions and check for last individual visit of a plant, within a time frame? (last one, for starters)


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
#                     # if (sp == orgs[o].sp && !(Base.isequal(org, orgs[o])) && org.reprd == false && orgs[o].reprd = false)
#                     #     #TODO add stochasticity
#                     #     org.reprd = true
#                     #     orgs[o].reprd = true
#                     #
#                     #     parents_genes = [org.genotype, orgs[o].reprd]
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
