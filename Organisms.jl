module Organisms

#using
using Distributions
using Setworld

#export
export newOrgs!, Organism, IDcounter

 #TODO check UNITS
const Boltz = 8.617e-5 # Boltzmann constant eV/K (non-SI) 1.38064852e-23 J/K if SI
const Ea = 0.69 # activation energy kJ/mol (non-SI), 0.63eV (MTE - Brown et al. 2004)
const Cgrowth = exp(25.2) # plant biomass production (Ernest et al. 2003) #TODO try a way of feeding those according to funcitonal group
const Cfertil = exp(26.0)
const Cmortal = exp(19.2)

"""
Organisms have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the rangge of search).
"""
mutable struct Organism
    ID::String
    #sp::String
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
    #Organism() = new()
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
    initorgs)
    # TODO optimize the for loops?
    orgs = []
    for frag in 1:size(landscape,3)
        for f in 1:length(initorgs.fgroups)
            # go through each x,y pair and place an org there
            #TODO might need outerconstructor for diffrent types of fgroups: insects dont have a radius, for example
            #for cell in eachindex(landscape[X,Y,frag])
            XYs = hcat(rand(1:size(landscape,1), initorgs.init_abund[f]),            rand(1:size(landscape,2), initorgs.init_abund[f])) #TODO indent it properly
            for i in 1:initorgs.init_abund[f]
                neworg = Organism( #TODO indent it tp Organsnim
                string(initorgs.fgroups[f], IDcounter + 1),
                [XYs[i,1] XYs[i,2] frag],
                initorgs.sps[f],
                initorgs.init_stage[f],
                0,
                false,
                initorgs.fgroups[f],
                ["placeholder" "placeholder"], #initialize with function
                rand(Distributions.Normal(initorgs.biomassμ[f],initorgs.biomasssd[f])),
                [initorgs.dispμ initorgs.dispshp],
                initorgs.radius[f])

                push!(orgs, neworg)
                #lanscape[X,Y,frag].orgs
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
    projmass!(landscape,orgs)
Projects the mass of each organisms stored in `orgs` into the `neighs` field of `landscape`. This projection means that the total biomass is divided into the square area delimited by the organism `radius`.

"""
function projvegmass!(landscape::Array{Any, N} where N,
    orgs::Array{Any,N} where N)

    for o in 1:length(orgs)
        x, y, frag = orgs[o].location
        r = orgs[o].radius
        fg = orgs[o].fgroup
        projmass = /(org[o].biomass, ((2*r+1)^2))

        for j in (y-r):(y+r), i in (x-r):(x+r)
            (i =< 0 || i > size(landscape[:,:,frag],1) || j =< 0 || j > size(landscape[:,:,frag],2)) && continue #check boundaries
            if haskey(landscape[i,j,frag].neighs,fg)
                landscape[i,j,frag].neighs[fg] += projmass # +1 because the center is a cell, not a point
            else
                landscape[i,j,frag].neighs[fg] = projmass
            end
        end
    end

end

# """
# """
function checkcompetition!(orgs::Array{Any,N} where N,
    landscape::Array{Any, N} where N)

    while o <= length(orgs)
        x, y, frag = orgs[o].location
        fg = orgs[o].fgroup
        #r = org.radius
        # 2. Look for neighbors in the area
        for j in (y-r):(y+r), i in (x-r):(x+r)
            (i == x && j == y) && continue #no self competition
            if landscape[i,j,frag].neighs[fg] > 0 # check the neighborhood of same fgroup for competition
                nbsum += landscape[i,j,frag].neighs
                nnbs += 1
            end
        end
        org[o].compterm = nbsum/nnbs # ratio of mass of neighs in nb of cells. TODO It should ne normalized somehow
        o += 1
    end
end

# """
# Individuals grow according to constraints from the MTE: local temperature affects biomass production rate. All individuals in a
# """
function grow!(orgs::Array{Any,N} where N,
    landscape,
    Ea,
    Boltz.,
    OrgsRef) #TODO organize this file to serve as refenrece of max size for different functional groups
    for o in 1:length(orgs)
        # dont grow beyond max size
        if orgs[o].biomass < #TODO check OrgsRef
            T = landscape[orgs[o].location].temp
            grown_mass += Cgrowth * organism.vegmass^(3/4) * exp(-Ea/(Boltz*T))
        else
            continue
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
