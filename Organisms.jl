module Organisms

"""
This module contains the

Organisms have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the rangge of search).
"""
#using
using Distributions
using Setworld

#export
export newOrgs!, Organism, IDcounter

 #TODO check UNITS
const Boltz = 8.617e-5 # Boltzmann constant eV/K (non-SI) 1.38064852e-23 J/K if SI
const aE = 0.69 # activation energy kJ/mol (non-SI), 0.63eV (MTE - Brown et al. 2004)
const Cgrowth = exp(25.2) # plant biomass production (Ernest et al. 2003) #TODO try a way of feeding those according to funcitonal group
const Cfertil = exp(26.0)
const tK = 273.15 # °C to K converter
const Cmortal = exp(19.2)

mutable struct Organism
    ID::String
    location::Array{Int64,2} # (fragment, x, y)#TODO use tuple
    sp::String
    stage::String #j,s,a #TODO seed = embryo
    age::Int64 # controls passing stages
    reped::Bool # reproduced?
    fgroup::String # Determines whether individual is a plant, insect, through functional group characterization
    genotype::Array{String,2}
    biomass::Float64
    disp::Array{Float64,2} #dispersal kernel parameters (mean and shape) TODO tuple
    radius::Int64
    Organism() = new()
end


"""
    newOrg(fgroups, init_abund, biomassμ, biomasssd)
newOrg() creates new `init_abund` individuals of each  functional group (`fgroups`), with biomass (`biomassμ` +- `biomasssd`). It is called during the initialization and all  the Organism's values are read from an parameter file for the especific functional group. All organisms, from all funcitonal groups, are stored in the same array.
#TODO those two take very different arguments. Still diffrentiate just methods?
    newOrg(parents, reprmode)
 newOrg() generates offspringn, after reproduction. For offspring creation, the `genotype` field is filled according to parental genotype and mode of reproduction and the other fields are copied from the parents.
# Methods
1. Upon initialization: `newOrg(fgroups, init_abund, biomassμ, landscape)`
`fgroups` list of funcitonal groups to be simulated
`init_abund` initial abundances
`biomassμ` mean biomass

# Arguments:
-`reprmode::String` specificies how the new organisms are generated:
    `init` initialization of the model
    `sex` sexual reproduction
    `clone` assexual reproduction
-`parent_s::Organism`
-`quant::Int64` is nb of new individuals or offspring to be created
"""

function newOrgs(landscape::Array{Any, N} where N,
                initorgs::InitOrgs)
    orgs = Organism[]
    for frag in 1:size(landscape,3)
        for f in 1:length(initorgs.fgroups) #TODO check the organisms file format: so far, all fragments get the same sps

            XYs = hcat(rand(1:size(landscape,1), initorgs.init_abund[f]),rand(1:size(landscape,2), initorgs.init_abund[f]))

            for i in 1:initorgs.init_abund[f]
                neworg = Organism(string(initorgs.fgroups[f], IDcounter + 1),
                                  [XYs[i,1] XYs[i,2] frag],
                                  initorgs.sps[f],
                                  initorgs.init_stage[f],
                                  0,
                                  false,
                                  initorgs.fgroups[f],
                                  ["placeholder" "placeholder"], #initialize with function
                                  rand(Distributions.Normal(initorgs.biomassμ[f],initorgs.biomasssd[f])),
                                  [initorgs.dispμ[f] initorgs.dispshp[f]],
                                  initorgs.radius[f])

                push!(orgs, neworg)

                global IDcounter += 1
            end
        end
    end
    #TODO different fgroups have different metabolic rates. Think if orgs should have a field for it  or what
    #TODO different stages have different biomass allocation for growht, reproduction, etc. Should there be a function to change those as stages change?
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
    orgs::Array{Any,N} where N)

    for o in 1:length(orgs)
        x, y, frag = orgs[o].location
        r = orgs[o].radius
        fg = orgs[o].fgroup
        projmass = /(org[o].biomass, ((2*r+1)^2))

        for j in (y-r):(y+r), i in (x-r):(x+r) #usar a funcao da FON Q trabalha com quadrantes?
            (i =< 0 || i > size(landscape[:,:,frag],1) || j =< 0 || j > size(landscape[:,:,frag],2)) && continue #check boundaries
            if fg == "autotroph" && i == x && j == y #"mark" steming point for plants
                landscape[i,j,frag].neighs[fg] = fgpars["autotrophsmax"]
            elseif haskey(landscape[i,j,frag].neighs,fg)
                landscape[i,j,frag].neighs[fg] += projmass
            else
                landscape[i,j,frag].neighs[fg] = projmass
            end
        end
    end

end

"""
    checkcompetition!()

"""
function checkcompetition!(orgs::Array{Any,N} where N, landscape::Array{Any, N} where N)

    while o <= length(orgs)
        x, y, frag = orgs[o].location
        fg = orgs[o].fgroup
        #r = org.radius
        # 2. Look for neighbors in the area
        for j in (y-r):(y+r), i in (x-r):(x+r)
            if landscape[i,j,frag].neighs[fg] > 0 # check the neighborhood of same fgroup for competition
                nbsum += landscape[i,j,frag].neighs
                nnbs += 1
            end
        end
        org[o].compterm = (nbsum - orgs[o].biomass)/nnbs # ratio of mass of neighs (- own biomass) in nb of cells. TODO It should ne normalized somehow
        o += 1
    end

end

"""
    allocate!()
Calculates biomass gain according to MTE rate and allocates it growth or reproduction, according to the developmental `stage` of the organism.
"""
function allocate!(orgs::Array{Any,N} where N,
    landscape,
    aE,
    Boltz.,
    OrgsRef) #TODO organize this file to serve as refenrece of max size for different functional groups

    for o in 1:length(orgs)

        T = landscape[orgs[o].location].temp + tK #(conversion to K)

        # Resource assimilation: #TODO check if resource allocation would be the same as growth
        # This MTE rate comes from dry weights: fat storage and whatever reproductive structures too, but not maintenance explicitly
        # Any cost related to insufficient minimal biomass goes into the survival probability function
        grown_mass += org[o].compterm * (Cgrowth * organism.vegmass^(3/4) * exp(-aE/(Boltz*T)))

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
    end

end

"""
    survive!()
Organism survival depends on total biomass, according to MTE rate.
"""
function survive!(orgs)
    indxs = []
    for o in 1:length(orgs)
        mort = Cmortal * (org[o].biomass["reprd"] + biomass["vegstruct"])^(-1/4)*exp(-aE/)
        if rand() > rand(PoissonBinomial(mort)) # successes are death events, because the consequence to modelled it to be taken out
            push!(indxs, o)
        end
    end
    deleteat!(orgs, indxs)
    return orgs
end

"""
    reproduce!()
Assigns proper reproduction mode according to the organism functional group.
"""
function reproduce!()
    #TODO expand choice of functional groups
    if "plant"
        reproduce!()
    elseif "insect"
        mate!()
    end

    offspring = []
    push!(offspring, Cfertil * org.biomass["reprd"]^(-1/4) * exp(-aE/(Boltz*T)))
end

"""
    mate!()
Insects reproduce if another one is found in the immediate vicinity.
"""
function mate!(org::Organism)
    x, y, frag = orgs.location #another org of same sp should match the locations of focus
    sp = org.sp

     # 1. check in the location field of orgs array:
     # 1.a inside same frag, look for locaions inside the squared area.
     # 2. when matching, differentiate between autotrphsa and the rest

    for o in 1:length(orgs) #look for partners
        #TODO optimize indexation of field location in arrray
        #TODO memory-wise, is it better to put all
        # 1:1 sex-ratio,
        if frag == orgs[o].location[3]
            if orgs[o].location[1:2] in collect(Iterators.product(x-1:x+1,y-1:y+1))
                #check sp, self and already reproduced
                if sp == orgs[o].sp && !(Base.isequal(org, orgs[o])) && org.reprd == false && orgs[o].reprd =
                    #TODO add stochasticity
                    org.reprd = true
                    orgs[o].reprd = true
                end
            end
        end
    end
end



# """
# Individuals reproduce with a fertility rate according to the MTE.
# When offspring is produced, they are generated
# """
# function reproduce!()
#
# end
#
#

end
