module Organisms

"""
This module contains the

Organisms have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the range of search).
"""
#using
using Distributions
using Setworld

#export
export newOrgs!, Organism, IDcounter

 #TODO CHECK units
const Boltz = 8.617e-5 # Boltzmann constant eV/K (non-SI) 1.38064852e-23 J/K if SI
const aE = 0.69 # activation energy kJ/mol (non-SI), 0.63eV (MTE - Brown et al. 2004)
const plants_gb0 = exp(25.2) # plant biomass production (Ernest et al. 2003) #TODO try a way of feeding those according to funcitonal group
const plants_fb0 = exp(26.0) # fertility rate
const tK = 273.15 # °C to K converter
const plants_mb0 = exp(19.2)

mutable struct Organism
    id::String
    location::Tuple{Int64} # (x,y,frag)
    sp::String
    stage::String #e,j,a
    age::Int64 # controls passing stages, phenology and
    reped::Bool # reproduced?
    fgroup::String # Plant, insect or more precise functional group
    genotype::Array{String,2}
    biomass::Dict
    disp::Tuple{Float64,2} #dispersal kernel parameters (a and b) TODO tuple
    radius::Array{Int64,2} # reproductive and vegetative area of influence. Not Tuple because not
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
`parent_s::Organism`
`quant::Int64` is nb of new individuals or offspring to be created
"""

function newOrgs(landscape::Array{Any, N} where N,
                initorgs::InitOrgs)
    orgs = Organism[]
    for frag in 1:size(landscape,3)
        for f in 1:length(initorgs.fgroups) #TODO check the organisms file format: so far, all fragments get the same sps

            XYs = hcat(rand(1:size(landscape,1), initorgs.init_abund[f]),rand(1:size(landscape,2), initorgs.init_abund[f]))

            for i in 1:initorgs.init_abund[f]
                neworg = Organism(string(initorgs.fgroups[f][1:3], IDcounter + 1),
                                  (XYs[i,1] XYs[i,2] frag),
                                  initorgs.sps[f],
                                  initorgs.init_stage[f],
                                  0,
                                  false,
                                  initorgs.fgroups[f],
                                  ["placeholder" "placeholder"], #initialize with function
                                  Dict("veg" => rand(Distributions.Normal(initorgs.biomassμ[f],initorgs.biomasssd[f]))),
                                  (initorgs.dispμ[f] initorgs.dispshp[f]),
                                  [initorgs.radius[f]])

                push!(orgs, neworg)

                global IDcounter += 1
            end
        end
    end
    #TODO different fgroups have different metabolic rates. Get it form OrgsRef
   return orgs #TODO or make sure that it has varied in the global scope as well rather export it?

end

# function newOrgs()
#     # TODO arguments are quite different. Should it be a different function?
#     # # TODO For the reproduction method of this fct: The genotype and individual location depends on the group of reproduction
# end

"""
    projvegmass!(landscape,orgs)
Projects the mass of each organisms stored in `orgs` into the `neighs` field of `landscape`. This projection means that the total biomass is divided into the square area delimited by the organism's `radius`.
"""
function projvegmass!(landscape::Array{Any, N} where N,
    orgs::Array{Organism,N} where N)

    for o in 1:length(orgs)
        x, y, frag = orgs[o].location
        r = orgs[o].radius
        fg = orgs[o].fgroup
        projmass = /(org[o].biomass["growth"], ((2*r+1)^2))

        for j in (y-r):(y+r), i in (x-r):(x+r) #usar a funcao da FON Q trabalha com quadrantes?
            if fg == "plant" && i == x && j == y #TODO center steming point has stronger biomass ~ way too close to FON again
                landscape[i,j,frag].neighs[fg] = fgpars["plant"]
            elseif haskey(landscape[i,j,frag].neighs,fg)
                landscape[i,j,frag].neighs[fg] += projmass
            else
                landscape[i,j,frag].neighs[fg] = projmass
            end
        end
    end

end

"""
    compete(orgs,landscape)
Each plant in `orgs` will check neighboring cells (inside it's zone of influence radius `r`) and detected overlaying ZOIs. When positive, the proportion of 'free' plant biomass is calculated: (focus plant biomass - sum(non-focus vegetative biomass))/(focus plant biomass), and normalized to the total area of projected biomass .
"""
function compete(org::Any,N where N, landscape::Array{Any, N} where N)

    x, y, frag = orgs[o].location
    fg = orgs[o].fgroup
    #r = org.radius
    # 2. Look for neighbors in the area
    for j in (y-r):(y+r), i in (x-r):(x+r) #TODO filter!() this area?
        if landscape[i,j,frag].neighs[fg] > 0 # check the neighborhood of same fgroup for competition
            nbsum += landscape[i,j,frag].neighs - orgs[o].biomass/((2*r+1)^2) #sum vegetative biomass of neighbors only (exclude focus plant own biomass)
        end
    end
    compterm = /(orgs[o].biomass - nbsum,orgs[o].biomass) # ratio of mass of neighs (- own biomass) in nb of cells. TODO It should ne normalized somehow

    return compterm
end



"""
    allocate!(orgs, landscape, aE, Boltz, OrgsRef)
Calculates biomass gain according to MTE rate and depending on competition. If competition is too strong, individual has a probability of dying. If not, gains is allocated to growth or reproduction, according to the developmental `stage` of the organism.
"""
function allocate!(orgs::Array{Any,N} where N,
    landscape,
    aE,
    Boltz,
    OrgsRef) #TODO organize this file to serve as refenrece of max size for different functional groups

    nogrowth = []

    for o in 1:length(orgs)

        T = landscape[orgs[o].location].temp + tK #(conversion to K)

        # Resource assimilation: #TODO check if resource allocation would be the same as growth
        # This MTE rate comes from dry weights: fat storage and whatever reproductive structures too, but not maintenance explicitly
        # Any cost related to insufficient minimal biomass goes into the survival probability function
        compterm = compete(orgs[o])
        if compterm > 0
            grown_mass += (1 - compterm) * (Cgrowth * organism.vegmass^(3/4) * exp(-aE/(Boltz*T)))

            #Resource allocation schedule
            #TODO make it more ellaborate and includde trade-offs
            if (org[o].stage == "e")  #TODO accountant for resistant stages?

                # embryos only consume reserves:
                org[o].biomass -=
            elseif stage == "j"
                # juveniles grow
                org[o].biomass["vegstruct"] += grown_mass
            elseif stage == "a"
                # adults reproduce
                org[o].biomass["reprd"] += grown_mass
            end
        end
    else
        push!(nogrowth,o)
    end
    return nogrowth
end

"""
    survive!(ors, nogrowth)
Organism survival depends on total biomass, according to MTE rate. However, the proportionality constants (b_0) used depend on the cause of mortality: competition-related, where
plants in nogrwth are subjected to two probability rates
"""
function survive!(orgs)

    deaths = []

    for o in 1:length(orgs)
        mortalrate = plants_mB #TODO call it from OrgsRef
        mB = mortalrate * (sum(values(orgs[o].biomass)))^(-1/4)*exp(-aE/)
        mprob = 1 - e^(-mB*orgs[o].age)

        # "weaker individuals" have two "ways" of dying
        if o in nogrowth
            compmortrate = plants_mB #TODO use different b_0 for mortality consequence of competition
            cmB = compmortrate * (sum(values(orgs[o].biomass)))^(-1/4)*exp(-aE/)
            cmprob =  1 - e^(-cmB*orgs[o].age)
            mprob += cmprob
        end

        if mprob > rand(PoissonBinomial(cmprob)) # successes are death events, because they are the value that is going to be relavant in here: the amount of individuals to be taken out
            push!(deaths, o)
        else
            orgs[o].age += 1
        end
    end
    deleteat!(orgs, deaths)
    return orgs
end

"""
    meanExP(a,b)
Draws a mean dispersed distance from the Exponential Power dispersal kernel.
"""
function meanExP(a::Float64,b::Float64)
    dist = Fileprep.lengthtocell(rand(a*(Distribution.Gamma(3/b))/(Distribution.Gamma(2/b))))
    return dist
end

"""
    checkboundaries(sourcefrag,dest)
`source` and `dest` contain the location indexes of the source (mother plant) and the pollen/seed. `checkboundaires()` verifies whether the new polen/seed location `(x,y)` is inside a habitat fragment (same as the source -`frag`- or another one insed the patch). Return a boolean that controls whether the process (reproduction or emergency/germination) proceeds or not.
"""
function checkboundaries(sourcefrag::Tuple{Int64}, xdest::Int64, ydest::Int64, fdest::Int64)
    #check inside frag

    if checkbounds(Bool, landscape[:,:,sourcefrag], xdest, ydest, fdest)
    #if (dx <= size(landscape[:,:,sf])[1] && dy <= size(landscape[:,:,sf])[2])
        inbound = true
    #elseif (dx <= size(landscape[:,:,sf])[1] && dy <= size(landscape[:,:,sf])[2])
    #TODO check ouside frag (more complex than checkbounds): how to detect the direction of neighboring fragments? use a matrix dependent probability?
        # -> compare dist with matriyes sizes and call possible fragmetns in the third index of location)
        #inbound = true
    else
        inbound = false
    end

    return inbound
end

"""
    reproduce!()
Assigns proper reproduction mode according to the organism functional group. This controls whether reproduction happens or not, for a given individual: plants depend on pollination, while insects do not. Following, it handles fertilization of new embryos and calculates offspring production. New individuals are included in the community at the end of current timestep.
"""
function reproduce!(orgs)
    #TODO sort out reproduction mode (pollination or not) according to funcitonal group
    # if # pollination depending plants
    #     pollination()
    # elseif "insect"
    #     parents_genes = mate!()
    # end

    reproducing <- filter.(x -> x.stage == "a", orgs)

    for o in 1:reproducing

        # Clonal reproduction

        # Find partners: Wind pollination
        θ = rand([0,0.5π,π,1.5π]) #get radian angle of distribution
        dist =  meanExP(2.3,0.44) #mean distance from the exponential power (Nathan et al. 2012), a = 2.3 and b = 0.44 according to Hardy et al. 2004.
        #new location
        xdest = round(Int64, reproducing[o] + dist*sin(theta), RoundNearestTiesAway)
        ydest = round(Int64, reproducing[o] + dist*cos(theta), RoundNearestTiesAway)
        fdest = reproducing[o].location[3] #TODO is landing inside the same fragment as the source, for now

        #check boundaries
        if checkboundaries(reproducing[o], xdest, ydest, fdest) #check if it lands somewhere

            # check for partners there #TODO make it less exact
            posptner = filter(x -> x.location .== (pollenland, reproducing[o]), orgs)

            if  length(posptner) => 1

                ptners = filter(x -> x.fgroup .== reproducing[o].fgroup)

                if length(ptners) > 1

                    #produce offsprings: newOrgs with it?
                    offspring = Organism[]
                    offsB = round(Cfertil * sum(values(reproducing[o].biomass))^(-1/4) * exp(-aE/(Boltz*T))) #TODO stochasticity!
                    for n in 1:noffsprg
                        #TODO check for a quicker way of creating several objects of composite-type
                        embryo = Organism(string(reproducing[o].fgroups, IDcounter + 1),
                        [], #location is given according to functional group and dispersal strategy, in disperse!()
                        reproducing[o].sp,
                        "e",
                        0,
                        false,
                        reproducing[o].fgroup,
                        ["placeholder" "placeholder"], #come from function
                        rand(Distributions.Normal(OrgsRef.biomassμ[reproducing[o].fgroup],OrgsRef.biomasssd[reproducing[o].fgroup])),
                        [OrgsRef.dispμ[f] OrgsRef.dispshp[f]],
                        OrgsRef.radius[f])
                        push!(offspring, embryo)
                    else
                        continue
                    end
                else
                    continue
                end
                # TODO pollination
            else
                continue
            end

        else
            continue

            # TODO gamete production and fertilization
            # parents_genes[1] and parents_genes[2]
        end

    end

end

"""
    disperse!()
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


"""
    mate!()
Insects reproduce if another one is found in the immediate vicinity.
"""
function mate!(org::Organism)
    x, y, frag = org.location #another org of same sp should match the locations of focus
    sp = org.sp

    # 1. check in the location field of orgs array:
    # 1.a inside same frag, look for locaions inside the squared area.
    # 2. when matching, differentiate between autotrphsa and the rest
    if org.stage == "a"

        for o in 1:length(orgs) #look for partners
            #TODO optimize indexation of field location in arrray
            #TODO memory-wise, is it better to put all ifs together?
            # 1:1 sex-ratio,
            if frag == orgs[o].location[3]
                if orgs[o].location[1:2] in collect(Iterators.product(x-1:x+1,y-1:y+1))
                    #check sp, self and already reproduced
                    if sp == orgs[o].sp && !(Base.isequal(org, orgs[o])) && org.reprd == false && orgs[o].reprd =
                        #TODO add stochasticity
                        org.reprd = true
                        orgs[o].reprd = true

                        parents_genes = [org.genotype, orgs[o].reprd]
                    end
                end
            end
        end

    else
        continue
    end
    return parents_genes #TODO check if it conflicts with modifying orgs
end


end
