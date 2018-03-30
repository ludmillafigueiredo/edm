"""
This module contains the

Organisms have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the range of search).
"""
module Organisms

using Distributions
using Setworld
using Fileprep

#export
export Organism, InitOrgs, newOrgs,  projvegmass!, allocate!, reproduce

 #TODO CHECK units
const Boltz = 8.617e-5 # Boltzmann constant eV/K (non-SI) 1.38064852e-23 J/K if SI
const aE = 0.69 # activation energy kJ/mol (non-SI), 0.63eV (MTE - Brown et al. 2004)
const plants_gb0 = exp(25.2) # plant biomass production (Ernest et al. 2003) #TODO try a way of feeding those according to funcitonal group
const plants_fb0 = exp(26.0) # fertility rate
const tK = 273.15 # °C to K converter
const plants_mb0 = exp(19.2)
const seedmassµ =

# Initial organisms parametrization
mutable struct InitOrgs
    fgroups::Tuple{String}
    sps::Tuple{String}
    init_stage::Tuple{String}
    init_abund::Tuple{Int64}
    #genotypes #TODO initialie those in functions
    biomassμ::Tuple{Float64}
    biomasssd::Tuple{Float64}
    dispμ::Tuple{Float64}
    dispshp::Tuple{Float64}
    radius::Tuple{Int64}
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
                                  (XYs[i,1], XYs[i,2], frag),
                                  initorgs.sps[f],
                                  initorgs.init_stage[f],
                                  0,
                                  false,
                                  initorgs.fgroups[f],
                                  ["placeholder" "placeholder"], #initialize with function
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
Projects the mass of each organisms stored in `orgs` into the `neighs` field of `landscape`. This projection means that the total biomass is divided into the square area delimited by the organism's `radius`.
"""
function projvegmass!(landscape::Array{Setworld.WorldCell, 3}, orgs::Array{Organism,N} where N)
    for o in 1:length(orgs)
        x, y, frag = orgs[o].location
        r = orgs[o].radius
        fg = orgs[o].fgroup
        projmass = /(orgs[o].biomass["veg"], ((2*r+1)^2))

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
function compete(landscape::Array{Setworld.WorldCell, 3},
    org::Organism)
    x, y, frag = org.location
    fg = org.fgroup
    r = org.radius
    # 2. Look for neighbors in the area
    nbsum = 0
    for j in (y-r):(y+r), i in (x-r):(x+r) #TODO filter!() this area?
        if !checkbounds(Bool,landscape[:,:,frag],j,i)
            continue
        else landscape[i,j,frag].neighs[fg] > 0 # check the neighborhood of same fgroup for competition
            nbsum += landscape[i,j,frag].neighs[fg] - org.biomass["veg"]/((2*r+1)^2) #sum vegetative biomass of neighbors only (exclude focus plant own biomass)
        end
    end
    compterm = /(org.biomass["veg"] - nbsum, org.biomass["veg"]) # ratio of mass of neighs (- own biomass) in nb of cells. TODO It should ne normalized somehow
    return compterm
end

"""
    allocate!(orgs, landscape, aE, Boltz, OrgsRef)
Calculates biomass gain according to MTE rate and depending on competition. If competition is too strong, individual has a higher probability of dying. If not, gains is allocated to growth or reproduction, according to the developmental `stage` of the organism.
"""
function allocate!(landscape::Array{Setworld.WorldCell,3},
    orgs::Array{Organism,N} where N,
    aE,
    Boltz) #TODO organize this object to serve as reference of max size for different functional groups

    nogrowth = Int64[]

    for o in 1:length(orgs)

        T = landscape[orgs[o].location[1], orgs[o].location[2], orgs[o].location[3]].temp

        # Resource assimilation: #TODO check if resource allocation would be the same as growth
        # This MTE rate comes from dry weights: fat storage and whatever reproductive structures too, but not maintenance explicitly
        # Any cost related to insufficient minimal biomass goes into the survival probability function
        compterm = compete(landscape, orgs[o])
        if compterm > 0
            grown_mass = (1 - compterm) * (plants_gb0 * sum(values(orgs[o].biomass))^(3/4) * exp(-aE/(Boltz*T)))

            #Resource allocation schedule
            #TODO make it more ellaborate and includde trade-offs
            if orgs[o].stage == "e"  #TODO accountant for resistant stages?
                # embryos only consume reserves: TODO realistic
                orgs[o].biomass["veg"] -= 0.01*org[o].biomass["veg"]
            elseif orgs[o].stage == "j"
                # juveniles grow
                orgs[o].biomass["veg"] += grown_mass
            elseif orgs[o].stage == "a" #TODO make it more complex. adults are investing everything in repoductiove biomass
                # adults reproduc
                if haskey(orgs[o].biomass,"reprd")
                    orgs[o].biomass["reprd"] += grown_mass
                else
                    orgs[o].biomass["reprd"] = grown_mass
                end
            end
        else
            push!(nogrowth,o)
        end
    end
    return nogrowth
end

"""
    meanExP(a,b)
Draws a mean dispersed distance from the Exponential Power dispersal kernel.
"""
function meanExP(a::Float64,b::Float64)
    dist = Fileprep.lengthtocell(a*rand(Distributions.Gamma(3/b))/rand(Distributions.Gamma(2/b)))
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
    reproduce(landscape,orgs)
Assigns proper reproduction mode according to the organism functional group. This controls whether reproduction happens or not, for a given individual: plants depend on pollination, while insects do not. Following, it handles fertilization of new embryos and calculates offspring production. New individuals are included in the community at the end of current timestep.
"""
function reproduce(landscape::Array{Setworld.WorldCell, 3}, orgs::Array{Organism,N} where N)
    #TODO sort out reproduction mode (pollination or not) according to funcitonal group
    # if # pollination depending plants
    #     pollination()
    # elseif "insect"
    #     parents_genes = mate!()
    # end

    reproducing = filter(x -> x.stage == "a", orgs)

    offspring = Organism[]

    for o in 1:length(reproducing)

        offsprgB = round(plants_fb0 * sum(values(reproducing[o].biomass))^(-1/4) * exp(-aE/(Boltz*T)),  RoundNearestTiesAway) #TODO stochasticity!

        # TODO there should be a realized number of offsprings that depends on the reproductive mass? offsprgB * seedmassµ

        for n in 1:offsprgB
            #TODO check for a quicker way of creating several objects of composite-type
            embryo = Organism(string(reproducing[o].fgroup[1:3], length(orgs) + length(offspring) + 1),
            (), #location is given according to functional group and dispersal strategy, in disperse!()
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

        reproducing[o].biomass["repr"] -= offsprgB*seedmass
    end

    return offspring
end

"""
    disperse!(offspring)
Butterflies, bees and seeds can/are disperse(d).
"""
function disperse!(orgs)
    # Seeds (Bullock et al. JEcol 2017)
    # Get seeds from org storage
    dispersing = filter.(x -> (x.fgroup = "plant" && x.stage == "e"), orgs) #TODO include insects
    # Exponential kernel for Ant pollinated herbs: mean = a.(Gamma(3/b)/Gamma(2/b))

    for d in dispersing
        # Dispersal distance from kernel:
        #acho q nao preciso de muito mais q Natahn et al. 2012 pra usar a pdf com as distancias da
        if d.fgroup == "wind"
            #Exp, herbs + appendage #TODO LogSech distribution
            dist = meanExP(4.7e-5,0.2336) #higher 9th percentile than ant (when comparing + appendage and 10-36 mg)
            # Check for available habitat (inside or inter fragment)
        elseif d.fgroup == "ant"
            #ExPherbs, 10-36 mg
            dist = meanExP(0.3726,1.1615)
        end
    end

    #TODO check border

end

"""
    survive!(ors, nogrowth,landscape)
Organism survival depends on total biomass, according to MTE rate. However, the proportionality constants (b_0) used depend on the cause of mortality: competition-related, where
plants in nogrwth are subjected to two probability rates
"""
function survive!(orgs::Array{Organism,N} where N, nogrowth::Array{Int64,N} where N,landscape::Array{Setworld.WorldCell,3})

    deaths = Int64[]

    for o in 1:length(orgs)

        T = landscape[orgs[o].location[1], orgs[o].location[2], orgs[o].location[3]].temp

        mortalconst = plants_mb0 #TODO call it from OrgsRef, when with different functional groups
        mB = mortalconst * (sum(values(orgs[o].biomass)))^(-1/4)*exp(-aE/T)
        mprob = 1 - e^(-mB*orgs[o].age)

        # individuals that didnt grow have
        if o in 1:length(nogrowth)
            compmortconst = plants_mb0 #TODO use different b_0 for mortality consequence of competition
            cmB = compmortconst * (sum(values(orgs[o].biomass)))^(-1/4)*exp(-aE/T)
            cmprob =  1 - e^(-cmB*orgs[o].age)
            mprob += cmprob
        end

        if rand() < mprob
            #mprob > rand(PoissonBinomial([mprob]),1)[1] # TODO should use a more ellaboratedistribution model? successes are death events, because they are the value that is going to be relavant in here: the amount of individuals to be taken out
            push!(deaths, o)
        else
            orgs[o].age += 1
        end
    end
    deleteat!(orgs, deaths)
    return orgs
end


"""
    emerge!()
Germination for seeds, emergency for eggs.
"""
function emerge!()
    embryos = filter.( x -> x.stage == "e", orgs)
    deaths = []

    for o in 1:length(embryos)
        if embryos[o].fgroup in ["plants"] #TODO sould be more dynamics and read from inputed fgroups
            germb0 = plants_gb0 #TODO call it from OrgsRe
        elseif embryos[o].fgroup in ["insects"]
            germb0 = plants_gb0 #TODO call it from OrgsRe
        end

        gB = germb0 * (sum(values(embryos[o].biomass)))^(-1/4)*exp(-aE/Boltz*T)
        gprob = 1 - e^(-gB*embryos[o].age) #TODO this age thing

        if gprob > rand(PoissonBinomial(gprob)) # successes are germinations
            embryos[o].stage = "j"
            embryos[o].age = 1
        else
            if (embryos[o].fgroup in ["plants"] &&  embryos[o].age <= 52) #TODO seed bank of one year. Make it more complex (count age)
                continue
            else
                push!(deaths, o)
            end
        end
    end
    deleteat!(orgs, deaths)
end

"""
    fitness(fgroup)
Calculates fitness for the functional group in
"""
function fitness()
end

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
