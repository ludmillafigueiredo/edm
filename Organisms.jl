module Organisms

#using
using Distributions
using Setworld

#export
export newOrgs!, Organism, IDcounter

const Boltz = 8.617e-5 # Boltzmann constant eV/K (non-SI) 1.38064852e-23 J/K if SI #TODO check SI
const Ea = 0.69 # activation energy kJ/mol (non-SI), 0.63eV (MTE - Brown et al. 2004)
const Cgrowth = exp(25.2) # plant biomass production (Ernest et al. 2003) #TODO try a way of feeding those according to funcitonal group
const Cfertil = exp(26.0)
const Cmortal = exp(19.2)
global IDcounter = Int64(0)

"""
Organisms have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the rangge of search).
"""
mutable struct Organism
    ID::String
    #sp::String
    location::Array{Int64,2} # (fragment, x, y)
    sp::String
    stage::String #j,s,a
    reped::Bool # reproduced?
    fgroup::String # Determines whether individual is a plant, insect, through functional group characterization
    genotype::Array{String,2}
    biomass::Float64
    disp::Array{Float64,2} #dispersal kernel parameters (mean and shape)
    radius::Int64
    #Organism() = new()
end


"""
    newOrg(fgroups, init_abund, biomassμ, biomasssd)
newOrg() creates new `init_abund` individuals of each  functional group (`fgroups`), with biomass (`biomassμ` +- `biomasssd`). It is called during the initialization and all  the Organism's values are read from an parameter file for the especific functional group.
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

function newOrgs!(landscape::Array{Setworld.WorldCell, N} where N,
    fgroups::Array{String, N} where N,
    sps::Array{String, N} where N,
    init_stage::Array{String, N} where N,
    init_abund::Array{Int64, N} where N,
    biomassμ::Array{Float64, N} where N,
    biomasssd::Array{Float64, N} where N,
    genotypes::Array{String, N} where N,
    dispμ::Array{Float64, N} where N,
    dispsd::Array{Float64, N} where N,
    radius::Array{Int64, N} where N,
    IDcounter::Int64)
    # TODO optimize the for loops
    orgs = []
    for frag in 1:size(landscape,3)
        for f in 1:length(fgroups)
            # go through each x,y pair and place an org there
            #TODO might need outerconstructor for diffrent types of fgroups: insects dont have a radius, for example
            #for cell in eachindex(landscape[X,Y,frag])
            XYs = hcat(rand(1:size(landscape,1), init_abund[f]),
            rand(1:size(landscape,2), init_abund[f]))
            for i in 1:init_abund[f]
                neworg = Organism(string(fgroups[f], IDcounter + 1),
                [XYs[i,1] XYs[i,2] frag],
                sps[f],
                init_stage[f],
                false,
                fgroups[f],
                [genotypes[f] "fillingERROR"],
                rand(srand(123), Distributions.Normal(biomassμ[f],biomasssd[f])),
                [dispμ dispsd],
                radius[f])

                push!(orgs, neworg)
                #lanscape[X,Y,frag].orgs
                IDcounter += 1
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

# """
# """
# function checkneighborhood(focalorga, landscape)
#     #pseudo:
#     #focalorga.location = x0,y0
#     #focalorga.radius = r
#     #neigh = landscape[x0-r:x0+r,y0-r:y0+r]
#     #check for organisms in that range: get it from landscape or from organims storage object?
# end

# """
# Individuals grow according to constraints from the MTE: local temperature affects biomass production rate. All individuals in a
# """
# function growth!(organism, landscape, Ea, Boltz)
#     grown_mass += Cgrowth * organism.vegmass^(3/4) * exp(-Ea/(Boltz*landscape[organism.location].temperature))
# end


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
