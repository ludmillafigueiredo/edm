
"""
This module contains the

Organisms have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the range of search).
"""
module Organisms

using Distributions
using JuliaDB
using DataValues
using StatsBase
using Setworld
using Fileprep

export Organism, OrgsRef, newOrgs!, projvegmass!, compete, develop!, allocate!, checkboundaries, reproduce!, mate!, mkoffspring!, disperse!, germinate, establish!, survive!, shedd!, destroyorgs!, release!

# Set up model constants
const Boltz = 8.62e-5 # eV/K Brown & Sibly MTE book chap 2
const aE = 0.65 # eV Brown & Sibly MTE book chap 2
const µ_wind = 0.1
const λ_wind = 3
const µ_ant = 1
const λ_ant = 0.2
#const plants_gb0 = (10^(10.15))/40 # 10e10.15 is the annual plant biomass production (Ernest et al. 2003) transformed to weekly base, with growth not happening during winter, and converted from kg to g
#const plants_mb0 = 9.902 #adjustted accordung to 1 death per individual for 1g for annuals (MTEpar notebook)
#const plants_fb0 = exp(30.0) # fertility rate

# Initial organisms parametrization
mutable struct OrgsRef
#species::Dict{String,String}
sp_id::Array{String, 1}
kernel::Dict{String,String}
kernel_sd::Dict{String,Float64}
e_mu::Dict{String,Float64}
e_sd::Dict{String,Float64}
b0g::Dict{String,Float64}
b0g_sd::Dict{String,Float64}
b0em::Dict{String,Float64}
b0em_sd::Dict{String,Float64}
b0am::Dict{String,Float64}
b0am_sd::Dict{String,Float64}
b0jg::Dict{String,Float64}
b0jg_sd::Dict{String,Float64}
b0ag::Dict{String,Float64}
b0ag_sd::Dict{String,Float64}
sestra::Dict{String,Int}
dyad::Dict{String,Float64}
floron::Dict{String,Int}
floron_sd::Dict{String,Int}
floroff::Dict{String,Int}
floroff_sd::Dict{String,Int}
seedon::Dict{String,Int}
seedon_sd::Dict{String,Int}
seedoff::Dict{String,Int}
seedoff_sd::Dict{String,Int}
max_mass::Dict{String,Float64}
first_flower::Dict{String,Int64}
first_flower_sd::Dict{String,Int64}
max_span::Dict{String,Int64}
max_span_sd::Dict{String,Int64}
mass_mu::Dict{String,Float64}
mass_sd::Dict{String,Float64}
abund::Dict{String,Int}
end

mutable struct Organism
id::String
location::Tuple # (x,y,frag)
sp::String #sp id, easier to read
mass::Dict
#### Evolvable traits ####
kernel::String
e_mu::Float64
b0g::Float64
b0em::Float64
b0am::Float64
b0jg::Float64
b0ag::Float64
#sestra::Int
#dyad::Float64
floron::Int
floroff::Int
seedon::Int
seedoff::Int
max_mass::Float64
first_flower::Int64
max_span::Int64
#mass_mu::Float64
#mass_sd::Float64
#### State variables ####
stage::String #e,j,a
age::Int64 # controls passing stages, phenology and
mated::Bool # TODO is it necessary?
genotype::Array{String,2} #initialize separately
radius::Int64
#Organism() = new()
end
Organism(id,location,sp,mass,kernel,e_mu,b0g,b0em,b0am,b0jg,b0ag,floron,floroff,seedon,seedoff,max_mass,first_flower,max_span) = Organism(id,location,sp,mass,kernel,e_mu,b0g,b0em,b0am,b0jg,b0ag,floron,floroff,seedon,seedoff,max_mass,first_flower,max_span,"a", 26,false,["A" "A"],0) #these individuals are initialized in the beginning of the simulation

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
function newOrgs!(landavail::Array{Bool,2},orgsref::Organisms.OrgsRef, id_counter::Int, tdist::String)

    orgs = Organism[]

    for frag in 1:size(landavail,3)

        for s in orgsref.sp_id

            XYs = hcat(rand(1:size(landavail,1),orgsref.abund[s]),
                       rand(1:size(landavail,2),orgsref.abund[s]))

            for i in 1:orgsref.abund[s]

                id_counter += 1 # update individual counter

                if tdist == "normal"
                    
                    neworg = Organism(hex(id_counter),
                                      (XYs[i,1],XYs[i,2],frag),
                                      s,
                                      Dict("veg" => 0.5*rand(Distributions.Normal(orgsref.mass_mu[s],orgsref.mass_sd[s])), #0.5 is aboveground biomass
                                           "repr" => 0),
                                      orgsref.kernel[s], #kernel
                                      rand(Distributions.Normal(orgsref.e_mu[s],orgsref.e_sd[s])), #e_mu
                    rand(Distributions.Normal(orgsref.b0g[s],orgsref.b0g_sd[s])), #b0g
                    rand(Distributions.Normal(orgsref.b0em[s],orgsref.b0em_sd[s])), #b0em
                    rand(Distributions.Normal(orgsref.b0am[s],orgsref.b0am_sd[s])), #b0am
                    rand(Distributions.Normal(orgsref.b0jg[s],orgsref.b0jg_sd[s])), #b0jg
                    rand(Distributions.Normal(orgsref.b0ag[s],orgsref.b0ag_sd[s])), #b0ag
                    Int(round(rand(Distributions.Normal(orgsref.floron[s],
                                                        orgsref.floron_sd[s])),RoundUp)), #floron
                    Int(round(rand(Distributions.Normal(orgsref.floroff[s],
                                                        orgsref.floroff_sd[s])),RoundUp)), #floroff
                    Int(round(rand(Distributions.Normal(orgsref.seedon[s],
                                                        orgsref.seedon_sd[s])),RoundUp)), #seedon
                    Int(round(rand(Distributions.Normal(orgsref.seedoff[s],
                                                        orgsref.seedoff_sd[s])),RoundUp)), #seedoff
                    Int(round(rand(Distributions.Normal(orgsref.max_span[s],
                                                        orgsref.max_span_sd[s])),RoundUp))) #max_span
                    
                    push!(orgs, neworg)                    
                else
                    neworg = Organism(hex(id_counter),
                                      (XYs[i,1],XYs[i,2],frag),
                                      s,
                                      Dict("veg" => 0.0,
                                           "repr" => 0.0),
                                      orgsref.kernel[s], #kernel
                                      rand(Distributions.Normal(mean(Uniform(min(orgsref.e_mu[s],orgsref.e_sd[s]),
                                                                             max(orgsref.e_mu[s],orgsref.e_sd[s])+0.00000001)),
                                                                abs(-(orgsref.e_mu[s],orgsref.e_sd[s])/6))), #e_mu
                    rand(Distributions.Normal(mean(Uniform(min(orgsref.b0g[s],orgsref.b0g_sd[s]),
                                                           max(orgsref.b0g[s],orgsref.b0g_sd[s])+0.00000001)),
                                              abs(-(orgsref.b0g[s],orgsref.b0g_sd[s])/6))), #b0g
                    rand(Distributions.Normal(mean(Uniform(min(orgsref.b0em[s],orgsref.b0em_sd[s]),
                                                           max(orgsref.b0em[s],orgsref.b0em_sd[s])+0.00000001)),
                                              abs(-(orgsref.b0em[s],orgsref.b0em_sd[s])/6))), #b0em
                    rand(Distributions.Normal(mean(Uniform(min(orgsref.b0am[s],orgsref.b0am_sd[s]),
                                                           max(orgsref.b0am[s],orgsref.b0am_sd[s])+0.00000001)),
                                              abs(-(orgsref.b0am[s],orgsref.b0am_sd[s])/6))), #b0am
                    rand(Distributions.Normal(mean(Uniform(min(orgsref.b0jg[s],orgsref.b0jg_sd[s]),
                                                           max(orgsref.b0jg[s],orgsref.b0jg_sd[s])+0.00000001)),
                                              abs(-(orgsref.b0jg[s],orgsref.b0jg_sd[s])/6))), #b0jg
                    rand(Distributions.Normal(mean(Uniform(min(orgsref.b0ag[s],orgsref.b0ag_sd[s]),
                                                           max(orgsref.b0ag[s],orgsref.b0ag_sd[s])+0.00000001)),
                                              abs(-(orgsref.b0ag[s],orgsref.b0ag_sd[s])/6))), #b0ag
                    Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.floron[s],orgsref.floron_sd[s]),
                                                                     max(orgsref.floron[s],orgsref.floron_sd[s])+0.00000001)),
                                                        abs(-(orgsref.floron[s], orgsref.floron_sd[s])/6))),RoundUp)), #floron
                    Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.floroff[s], orgsref.floroff_sd[s]),
                                                                     max(orgsref.floroff[s], orgsref.floroff_sd[s])+0.00000001)),
                                                        abs(-(orgsref.floroff[s], orgsref.floroff_sd[s])/6))),RoundUp)), #floroff
Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.seedon[s], orgsref.seedon_sd[s]),
                                                 max(orgsref.seedon[s], orgsref.seedon_sd[s])+0.00000001)),
                                    abs(-(orgsref.seedon[s],orgsref.seedon_sd[s])/6))),RoundUp)), #seedon
Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.seedoff[s], orgsref.seedoff_sd[s]),
                                                 max(orgsref.seedoff[s], orgsref.seedoff_sd[s])+0.00000001)),
                                    abs(-(orgsref.seedoff[s],orgsref.seedoff_sd[s])/6))),RoundUp)), #seedoff
0.0,
Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.first_flower[s], orgsref.first_flower_sd[s]),
                                                 max(orgsref.first_flower[s], orgsref.first_flower_sd[s])+0.00000001)),
                                    abs(-(orgsref.max_span[s], orgsref.first_flower[s])/6))),RoundUp)),
Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.max_span[s], orgsref.max_span_sd[s]),
                                                 max(orgsref.max_span[s], orgsref.max_span_sd[s])+0.00000001)),
                                    abs(-(orgsref.max_span[s], orgsref.max_span_sd[s])/6))),RoundUp)))#max_span

neworg.max_mass = (neworg.e_mu*1000/2.14)^2
neworg.mass["veg"] = neworg.max_mass * 0.5

push!(orgs, neworg)
end

end
end
end

return orgs, id_counter

end

"""
projvegmass!(landscape,orgs)
Rewrites the projected mass of each organisms stored in `orgs` into the `neighs` field of `landscape`. This projection means that the total biomass is divided into the square area delimited by the organism's `radius`.
"""
function projvegmass!(landscape::Array{Dict{String, Float64},2}, orgs::Array{Organism,1}, settings::Dict{String, Any})
    
    competing = find(x->(x.stage == "a" || x.stage == "j"),orgs) #juveniles com ashard as adults, but have higher growth rate and lower mortality

    for o in competing
        x, y, frag = orgs[o].location
        orgs[o].radius = round(Int64, (1/4) * (sqrt(orgs[o].mass["veg"]^(2/3)) - 1)/2, RoundUp) # multiply by 0.2 because weiner uses 1cm2 projections, and my cells are 16 cm2

        if orgs[o].radius == 0  #TODO: hotfix, check a more sound solution: vegetative mass - 1 might be too much for juveniles, especially the young ones
            orgs[o].radius = 1
        end
        
        r = orgs[o].radius # separated for debugging

        projmass = /(orgs[o].mass["veg"], ((2*r+1)^2))

        # unitytest
        open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            println(sim, orgs[o].id," has ",orgs[o].mass["veg"], " radius $r and projvegmass:", projmass) #ugly format to avoid risking some anoying errors that have been happening
        end

        for j in (y-r):(y+r), i in (x-r):(x+r) #TODO usar a funcao da FON Q trabalha com quadrantes? dar mais peso para steming point?
            if !checkbounds(Bool,landscape[:,:,frag],j,i) # check boundaries: absorbing borders: the biomass is not re-divided to the amount of cells inside the fragment. What is projected outside the fragmetn is actually lost: Edge effect
                continue
            else
                if haskey(landscape[i,j,frag],"p")
                    
                    landscape[i,j,frag]["p"] += projmass
                    # unity test
                    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                        println(sim, orgs[o].id,"is in $(orgs[o].location) and  has found a neighbor.")
                        println(sim, "Already projected in x = $i, y = $j:", landscape[i,j,frag]["p"])
                    end
                    
                else
                    landscape[i,j,frag] = Dict("p" => projmass)
                    # unity test
                    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                        println(sim, orgs[o].id," has no neighbor")
                    end
                    
                end
            end
        end
    end
end

"""
compete(landscape, orgs)
Each plant in `orgs` will check neighboring cells (inside it's zone of influence radius `r`) and detected overlaying ZOIs. When positive, the proportion of 'free' plant biomass is calculated: (focus plant biomass - sum(non-focus vegetative biomass))/(focus plant biomass), and normalized to the total area of projected biomass .
"""
function compete(landscape::Array{Dict{String, Float64},2}, org::Organism,settings::Dict{String, Any})

    x, y, frag = org.location
    sp = org.sp
    r = org.radius

    compterm = 0
    nbsum = 0

    # 2. Look for neighbors in the square area (2r+1) delimited by r
    for j in (y-r):(y+r), i in (x-r):(x+r)
        if !checkbounds(Bool,landscape[:,:,frag],i,j)
            continue
            # #unity test
            open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                println(sim, "$(org.id) out of bound projection at $(org.location)")
                #println("$(org.id) out of bound projection at $(org.location)")
            end
        elseif j == y && i == x # steming point is "stronger", doesnt compete
            continue
        elseif haskey(landscape[i,j,frag],"p")
            #landscape[i,j,frag].neighs[fg] > 0 # check the neighborhood of same fgroup for competition
            nbsum += (landscape[i,j,frag]["p"] - /(org.mass["veg"],(2*r+1)^2)) #sum vegetative biomass of neighbors only (exclude focus plant own biomass projection in that cell)
            # unity test
            open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                #println(org.id," has nbsum = ", nbsum, "from $i, $j")
                println(sim, org.id," has nbsum = ", nbsum) #ugly format to avoid risking some anoying errors that have been happening
            end
        end
    end

    compterm = /((org.mass["veg"] - nbsum), org.mass["veg"])
    # unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        #println("$(org.id) radius = $r and compterm = $compterm")
        println(sim, "$(org.id) radius = $r and compterm = $compterm")
    end
    
    return compterm
end


"""
allocate!(orgs, landscape, aE, Boltz, OrgsRef)
Calculates biomass gain according to MTE rate and depending on competition. Competition is measured via a biomass-based index `compterm` (`compete` function). This term gives the proportion of actual biomass gain an individual has. If competition is too strong (`compterm` < 0), the individual has a higher probability of dying. If not, the biomass is allocated to growth or reproduction, according to the developmental `stage` of the organism and the season (`t`) (plants start allocating to week 12).
"""
function allocate!(landscape::Array{Dict{String, Float64},2}, orgs::Array{Organism,1}, t::Int64, aE::Float64, Boltz::Float64, settings::Dict{String, Any},orgsref::Organisms.OrgsRef,T::Float64)
    #1. Initialize storage of those that are ont growing and have higher prob of dying (later)
    nogrowth = Int64[]

    if settings["competition"] != "capacity"
        growing = find(x -> x.stage in ("a","j"), orgs)
        for o in growing

            # This MTE rate comes from dry weights: fat storage and whatever reproductive structures too, but not maintenance explicitly
            # Any cost related to insufficient minimal biomass goes into the survival probability function
            # 2.b Check for competition

            compterm = compete(landscape, orgs[o], settings)
            # unity test
            #println(simulog, org.id," weights",org.biomass["veg"]," had $nbsum g overlap")
            #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #    println(sim, orgs[o].id, "  compterm $compterm")
            #end
            
            # ZOI based competition
            if compterm <= 0
                #2.c Those not growing will have higher chance of dying
                push!(nogrowth,o)
            else
                
            end
        end
    else
        # seedling and adults grow at different rates
        competing = find(x->(x.stage == "a" || x.stage == "j"),orgs)
        for o in competing

            if orgs[o].stage == "j"
                b0 = orgs[o].b0jg
            elseif orgs[o].stage == "a"
                b0 = orgs[o].b0ag
            else
                error("Seed or egg trying to compete.")
            end

            #println("$(sum(collect(values(orgs[o].mass))))")
            grown_mass = b0*(sum(collect(values(orgs[o].mass))))^(3/4)*exp(-aE/(Boltz*T))

            if grown_mass <= 0
                push!(nogrowth,o)
            end            
            
            # unity test
            #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #    println(sim, "$(orgs[o].id) should grow $grown_mass")
            #end

            #Resource allocation schedule
            #TODO make it more ellaborate and includde trade-offs
            
            if orgs[o].stage == "j"
                # juveniles grow
                orgs[o].mass["veg"] += grown_mass 
                # unity test
                #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a")do sim
                #    println(sim, "$(orgs[o].id)-$(orgs[o].stage) grew $grown_mass")
                #end
            elseif orgs[o].stage == "a" &&
                (orgs[o].floron <= rem(t,52) < orgs[o].floroff) && (sum(collect(values(orgs[o].mass))) >= 0.1*(orgs[o].max_mass))
                # adults invest in reproduction
                if haskey(orgs[o].mass,"repr")
                    orgs[o].mass["repr"] += grown_mass 
                    #unity test
                    #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                    #    println(sim, "$(orgs[o].id)-$(orgs[o].stage) is FLOWERING repr= ",orgs[o].mass["repr"])
                    #end
                else
                    orgs[o].mass["repr"] = grown_mass
                    #unity test
                    #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                    #    println(sim, "$(orgs[o].id)-$(orgs[o].stage) is FLOWERING repr= ",orgs[o].mass["repr"])
                    #end
                end
            elseif orgs[o].stage == "a" && sum(values(orgs[o].mass)) < orgs[o].max_mass # TODO refer it to a 50% of the species biomass #individuals that are too small dont reproduce
                orgs[o].mass["veg"] += grown_mass 
                # unity test
                #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                #    println(sim, "$(orgs[o].id)-$(orgs[o].stage) grew VEG $grown_mass")
                #end
                #else
                #continue #error("Seed or egg trying to allocate.")
            end

            #unity test
            # open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #     println(sim, "current biomass: $(orgs[o].biomass)")
            # end
            
        end
  
        # unity test
        #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        #    println(sim, "Not growing $nogrowth")
        #end
    end
end
"""
develop!()
Controls individual stage transition.
"""
function develop!(orgs::Array{Organism,1}, orgsref::Organisms.OrgsRef)
    juvs = find(x->x.stage == "j",orgs)

    for j in juvs
        if orgs[j].age >= orgs[j].first_flower
            orgs[j].stage = "a"
        end
    end
    #TODO insect holometabolic growth
    #TODO insect hemimetabolic development
end

"""
checkboundaries(sourcefrag,xdest, ydest, fdest)
`source` and `dest` contain the location indexes of the source (mother plant) and the pollen/seed. `checkboundaires()` verifies whether the new polen/seed location `(x,y)` is inside a habitat fragment (same as the source -`frag`- or another one insed the patch). Return a boolean that controls whether the process (reproduction or emergency/germination) proceeds or not.
"""
function checkboundaries(landavail::Array{Bool,2}, xdest::Int64, ydest::Int64, fdest::Int64) #TODO is this function necessary?
    #check inside frag
    if checkbounds(Bool, landavail[:,:,fdest], xdest, ydest)
        inbound = true
    else
        # if LandPars.nfrags > 1
        #  if dist in connectivity matrix
        #   xdest = rand(x from connectmatrix[curr,dest])
        #   ydest = rand(y from connect_matrix[current,dest])
        #   fdest = connect column
        #   inbound = true
        #  end
        # else
        inbound = false
        # end
    end
    return inbound #,xdest, ydest, fdest
end

"""
reproduce!(landscape,orgs)
Assigns proper reproduction mode according to the organism functional group. This controls whether reproduction happens or not, for a given individual: plants depend on pollination, while insects do not. Following, it handles fertilization of new embryos and calculates offspring production. New individuals are included in the community at the end of current timestep.
    Will probably call mate!()
    """
function reproduce!()
    
end

"""
    mate!()
    Calculate proportion of insects that reproduced (encounter?) and mark that proportion of the population with the `mated` label.
    GEnetic comes here
    """
function mate!(orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any}, scen, tp, remain, regime, kp)

    ready = find(x-> x.stage == "a" && x.mass["repr"] > x.e_mu, orgs) # TODO find those with higher reproductive mas than the mean nb of seeds * seed mass.
    pollinated = []

    # POLLINATION INDEPENDENT SCENARIO: as long as there is reproductive allocation, there is reproduction. Clonality is virtually useless
    if length(ready) > 0 # check if there is anyone flowering 
        
        if scen == "indep"
            for r in ready
                orgs[r].mated = true
            end

            # unity test
            #open(string("EDoutputs/",settings["simID"],"simulog.txt"), "a") do sim
            #println("Pollination scenario: $scen")
            #end
            
            # POLLINATION DEPENDENT SCENARIO:
        elseif scen == "equal"

            # unity test
            #open(string("EDoutputs/",settings["simID"],"simulog.txt"), "a") do sim
            #println("Pollination scenario: $scen")
            #end

            if t < tp
                pollinated = ready
                for p in pollinated 
                    orgs[p].mated = true
                end
            else 
                if regime == "pulse"
                    if t == tp # pollination lost only once
                        pollinated = sample(ready,n=Int(floor(length(ready)* exp(-(t-tp)) * kp)), replace = false, ordered = true) #* 10^(-3)) 
                    else
                        pollinated = ready
                    end

                    for p in pollinated
                        orgs[p].mated = true
                    end
                    
                elseif regime == "exp"
                    pollinated = sample(ready, Int(floor(length(ready)* exp(-(t-tp)) * kp)), replace = false, ordered = true) #* 10^(-3)) # (kp = 1) rdmly take the nbs of occupied flowers/individuals to be set to reproduce#
                    # exp(tp - t) makes the pollination loss decrease from 1 (tp = t) to 0
                    for p in pollinated
                        orgs[p].mated = true
                    end
                elseif regime == "const"
                    pollinated = sample(ready, Int(floor(length(ready)* remain * kp)), replace = false, ordered = true) #* 10^(-3))
                    # unity test
                    open(string("EDoutputs/",settings["simID"],"simulog.txt"), "a") do sim
                        println(sim,"Pollinated in $regime regime: $pollinated")
                    end 
                    for p in pollinated
                        orgs[p].mated = true
                        # unity test
                        open(string("EDoutputs/",settings["simID"],"simulog.txt"), "a") do sim
                            println(sim,"$(orgs[p].id) can reproduce? $(orgs[p].mated)")
                        end 
                    end
                else
                    error("Please chose a pollination scenario \"scen\" in insect.csv:
                          - \"pulse\":
                          - \"const\": or
                          - \"exp\"")                    
                end
                
            end
            
            # elif settings["insect"] == "spec"
            #   ready = rand(length(repro)* remain * kp * 10-³), repro)
        else
            error("Please chose a pollination scenario \"scen\" in insect.csv:
                  - \"indep\": pollination independent scenario,
                  - \"depend\": pollination dependent scenario, with supplementary info for simulating.")
        end
    end
    
    #unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println(sim, "Reproducing: $(length(ready)), week $t.
                Actually reprod: $pollinated , because $scen and $regime.")
    end
    
end

"""
    mkoffspring!()
After mating happened (marked in `reped`), calculate the amount of offspring
 """
function mkoffspring!(orgs::Array{Organisms.Organism,1}, t::Int64, settings::Dict{String, Any},orgsref::Organisms.OrgsRef, id_counter::Int)
    
    offspring = Organism[]

    for sp in unique(getfield.(orgs, :sp))

        ferts = find(x -> x.mated == true && x.sp == sp, orgs)
        
        #unity test
        #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        #println("Reproducing: $ferts week $t")
        #end

        for o in ferts

            emu = orgs[o].e_mu
            offs =  div(orgs[o].mass["repr"], emu)
            # unity test
            open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                println(sim, orgs[o].id, " has biomass $(orgs[o].mass["repr"]) and produces $offs seeds.")
            end

            if offs <= 0
                continue
            else
                
                orgs[o].mass["repr"] -= (offs * emu)

                # unity test
                open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                    println(sim, orgs[o].id, " now has biomass $(orgs[o].mass["repr"]).")
                end

                if  orgs[o].mass["repr"] <= 0
                    orgs[o].mass["repr"] = 0
                end

                #unity test
                #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                #println("Number of offspring of ", orgs[o].id, ": ",offs)
                #end

                # get phylogenetic constraint: variance of the distribution
                sp = orgs[o].sp
                conspp = [rand(filter(x -> x.sp == sp && x.stage == "a", orgs))]

                e_mudist = Array{Float64}(1,length(conspp)) 
                b0gdist = Array{Float64}(1,length(conspp))
                b0emdist = Array{Float64}(1,length(conspp))
                b0amdist = Array{Float64}(1,length(conspp))
                b0jgdist = Array{Float64}(1,length(conspp))
                b0agdist = Array{Float64}(1,length(conspp))
                florondist = Array{Int}(1,length(conspp))
                floroffdist = Array{Int}(1,length(conspp))
                seedondist = Array{Int}(1,length(conspp))
                seedoffdist = Array{Int}(1,length(conspp))
                first_flowerdist = Array{Float64}(1,length(conspp))
                max_spandist = Array{Int64}(1,length(conspp))

                for i in eachindex(conspp)
                    e_mudist[i] = conspp[i].e_mu 
                    b0gdist[i] = conspp[i].b0g
                    b0emdist[i] = conspp[i].b0em
                    b0amdist[i] = conspp[i].b0am
                    b0jgdist[i] = conspp[i].b0jg
                    b0agdist[i] = conspp[i].b0ag
                    florondist[i] = conspp[i].floron
                    floroffdist[i] = conspp[i].floroff
                    seedondist[i] = conspp[i].seedon
                    seedoffdist[i] = conspp[i].seedoff
                    first_flowerdist[i] = conspp[i].first_flower
                    max_spandist[i] = conspp[i].max_span
                end
                
                for n in 1:offs
                    
                    id_counter += 1 # update individual counter

                    embryo = deepcopy(orgs[o])

                    newvalue = rand(Distributions.Normal(0,abs(embryo.e_mu-mean(e_mudist) + 0.00000001)))
                    embryo.e_mu + newvalue >= orgsref.e_mu[embryo.sp] ? # if seed biomass or minimal biomass would smaller than zero, it does not chenge
                    embryo.e_mu += newvalue : embryo.e_mu += 0
                    embryo.b0g += rand(Distributions.Normal(0,abs(embryo.b0g-mean(b0gdist) + 0.00000001)/embryo.b0g))
                    embryo.b0em += rand(Distributions.Normal(0,abs(embryo.b0em-mean(b0emdist) + 0.00000001)/embryo.b0em))
                    embryo.b0am += rand(Distributions.Normal(0,abs(embryo.b0am-mean(b0amdist) + 0.00000001)/embryo.b0am))
                    embryo.b0jg += rand(Distributions.Normal(0,abs(embryo.b0jg-mean(b0jgdist) + 0.00000001)/embryo.b0jg))
                    embryo.b0ag += rand(Distributions.Normal(0,abs(embryo.b0ag-mean(b0agdist) + 0.00000001)/embryo.b0ag))
                    embryo.floron += Int(round(rand(Distributions.Normal(0,abs(embryo.floron-mean(florondist) + 0.00000001)/embryo.floron)),RoundUp))
                    embryo.floroff += Int(round(rand(Distributions.Normal(0,abs(embryo.floroff-mean(floroffdist) + 0.00000001)/embryo.floroff)),RoundUp))
                    embryo.seedon += Int(round(rand(Distributions.Normal(0,abs(embryo.seedon-mean(seedondist) + 0.00000001)/embryo.seedon)),RoundUp))
                    embryo.seedoff += Int(round(rand(Distributions.Normal(0,abs(embryo.seedoff-mean(seedoffdist) + 0.00000001)/embryo.seedoff)),RoundUp))
                    newvalue = Int(round(rand(Distributions.Normal(0,abs(embryo.first_flower-mean(first_flowerdist) + 0.00000001)/embryo.first_flower)), RoundUp))
                    embryo.first_flower + newvalue > 12 ?
                    embryo.first_flower += newvalue : embryo.first_flower += 0
                    newvalue = Int(round(rand(Distributions.Normal(0,abs(embryo.max_span-mean(max_spandist) + 0.00000001)/embryo.max_span)),RoundUp))
                    embryo.max_span + newvalue > orgsref.max_span[embryo.sp] ?
                    embryo.max_span += newvalue : embryo.max_span += 0

                    # reset adult max_mass according to newly set 
                    embryo.max_mass = (embryo.e_mu*1000/2.14)^2
                    
                    # set embryos own individual non-evolutionary traits
                    embryo.id = hex(id_counter)
                    embryo.location = orgs[o].location #stays with mom until release 
                    embryo.mass = Dict("veg" => embryo.e_mu,
                                       "repr" => 0.0)
                    embryo.stage = "e"
                    embryo.age = 0
                    embryo.mated = false
                    embryo.genotype = ["A" "A"]
                    embryo.radius = 0

                    push!(offspring, embryo)
	            #unity test
                    #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a")do sim
                    #println("Pushed new org ", embryo, " into offspring")
                    #end
                end
orgs[o].mated = false # after producing seeds in a week, the plant will only do it again in the next week if it gets pollinated again
end
end

end

append!(orgs, offspring)

#unity test
open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
    println(sim, "Total offspring at week $t: " ,length(offspring))
end

return id_counter

end

"""
release!()
Probably need to call disperse here, because not all "e"s in orgs are released at the same time.
"""
function release!(orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any},orgsref::Organisms.OrgsRef )
    # Individuals being released in any given week are: in embryo stage (=seed= & in their seed release period (seedon <= t <= seedoff for the species)
    seedsi = find(x -> x.stage == "e" && x.age == 0 && x.seedon <= rem(t,52) < x.seedoff, orgs) # using a condition "outside" orgs might not work. This condition with orgsref only works because orgsref always has the sp names of x as keys in the dictionnary. If presented with a key that it does do contain, it throws an error.
    return seedsi
end

"""
disperse!(landscape, orgs, t, seetings, orgsref,)
Seeds are dispersed.
"""
function disperse!(landavail::Array{Bool,2},seedsi, orgs::Array{Organisms.Organism, 1}, t::Int, settings::Dict{String, Any},orgsref::Organisms.OrgsRef, connects::Array{Float64,2}, AT::Float64, Ah:Float64)

    #unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println(sim, "Dispersing: $seedsi")
    end

    lost = Int64[]

    for d in seedsi
        if orgs[d].kernel == "a"
            dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ_ant,λ_ant),1)[1])
            #unity test
            #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #    println(sim, "$(orgs[d].id) ant dispersal kernel")
            #end
        elseif orgs[d].kernel == "w"
            dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ_wind,λ_wind),1)[1])
            #unity test
            #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #    println(sim, "$(orgs[d].id) wind dispersal kernel")
            #end
        elseif orgs[d].kernel == "wa"
            µ,λ = rand([[µ_wind λ_wind],[µ_ant λ_ant]])
            dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ,λ),1)[1])
            #unity test
            #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #    println(sim, "$(orgs[d].id) $µ dispersal kernel")
            #end
        else
            error("Check dispersal kernel input for species $(orgs[d].sp).")
        end
        
        #unity test
        #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println("$(orgs[d].id) dispersal distance: $dist")
        #end

        # Find patch
        θ = rand([0 0.5π π 1.5π]) # TODO tentar rand(0:2)*pi?
        xdest = orgs[d].location[1] + dist*round(Int64, cos(θ), RoundNearestTiesAway)
        ydest = orgs[d].location[2] + dist*round(Int64, sin(θ), RoundNearestTiesAway)
        fsource = orgs[d].location[3] 

        if checkboundaries(landavail, xdest, ydest, fsource) && landavail[xdest, ydest, fsource] == true # intra fragment dispersal
            orgs[d].location = (xdest,ydest,fsource)
            #unity test
            #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #println(orgs[d].id," dispersed $dist")
            #end
        else # out of fragment. Try another one
            # check if someone is inside the distance:
            fdest = find(x -> x <= dist && x > 0, connects[:,fsource]) |> (y -> (length(y) > 0 ? fdest = rand(x) : false))
            if fdest == false
                push!(lost,d)
            else
                rand(Distributions.Binomial(1,AT/Ah),1)[1] == 1 ? orgs[d].location = (xdest,ydest,fdest) : push!(lost,d)
            end
            #unity test
            #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #println(orgs[d].id," dispersed $dist but died")
            #end
        end
        # else
        #unity test
        # open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        #     println(sim, orgs[d].id, orgs[d].sp," Not releasing seeds.")
        # end
        # end
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
function germinate(org::Organisms.Organism)
    gprob = 1 - exp(-org.b0jg)
    if gprob < 0
        gprob = 0
    elseif gprob > 1
        gprob = 1
    end
    
    germ = false
    if 1 == rand(Distributions.Binomial(1,gprob),1)[1]
        germ = true
    end
    return germ
end

"""
establish!
Seeds only have a chance of establishing in patches not already occupied by the same funcitonal group, in. When they land in such place, they have a chance of germinating (become seedlings - `j` - simulated by `germinate!`). Seeds that don't germinate stay in the seedbank, while the ones that are older than one year are eliminated.
    """
function establish!(landscape::Array{Dict{String,Float64},2}, orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any}, orgsref::Organisms.OrgsRef)
    #REFERENCE: May et al. 2009
    establishing = find(x -> x.stage == "e", orgs)

    #unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println(sim, "$(length(establishing)) seed trying to establish")
    end

    lost = Int64[]

    for o in establishing
        if orgs[o].seedon <= rem(t,52) #only after release and dispersal a seed can establish
            orgcell = orgs[o].location
            if haskey(landscape[orgcell[1], orgcell[2], orgcell[3]],["p"])
                push!(lost,o)
	        #unity test
   	        open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        	    println(sim, "$o falling into occupied cell and died")
    	        end
	    end

            if germinate(orgs[o])
                orgs[o].stage = "j"
	        #unity test
   	        open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        	    println(sim, "Became juvenile: $(orgs[o].id)")
    	        end
            end
        end
        
    end
    deleteat!(orgs,lost)
end

"""
    survive!(orgs, nogrowth,landscape)
    Organism survival depends on total biomass, according to MTE rate. However, the proportionality constants (b_0) used depend on the cause of mortality: competition-related, where
    plants in nogrwth are subjected to two probability rates
    """
function survive!(orgs::Array{Organisms.Organism,1}, t::Int, currentk::Float64, settings::Dict{String, Any}, orgsref::Organisms.OrgsRef, landavail::Array{Bool,2},T)

    deaths = Int64[]
    seeds = find(x -> x.stage == "e", orgs) 

    # Density-independent mortality
    for o in 1:length(orgs)

        #T = landscape[orgs[o].location[1], orgs[o].location[2], orgs[o].location[3]].temp

        if sum(values(orgs[o].mass)) <= 0
            mprob = 1
        elseif o in seeds
            if rem(t,52) > orgs[o].seedoff #seeds that are still in the mother plant cant die. If their release season is over, it is certain thatthey are not anymore, even if they have not germinated 
                Bm = orgs[o].b0em * (orgs[o].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
                println("Bm: $Bm, b0em = $(orgs[o].b0em), seed mass = $(orgs[o].mass["veg"])")
                mprob = 1 - exp(-0.1)
            else
                continue
            end 
        elseif (orgs[o].age/52) >= orgs[o].max_span #oldies die
            mprob = 1
            #unity test
            open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                println(sim, "$(orgs[o].id) $(orgs[o].stage) dying of old age")
            end
        else
            Bm = orgs[o].b0am * (sum(collect(values(orgs[o].mass))))^(-1/4)*exp(-aE/(Boltz*T))
            # unity test
            #open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
             #   println(sim,"$(orgs[o].id) $(orgs[o].stage) mortality rate $Bm")
            #end
            mprob = 1 - exp(-Bm)               
        else
            mprob = 0            
        end

        if mprob < 0
            mprob = 0
        elseif mprob > 1
            mprob = 1
        end
        
        if 1 == rand(Distributions.Binomial(1,mprob),1)[1] || sum(values(orgs[o].mass)) <= 0
            currentk -= sum(values(orgs[o].mass)) #update currentk to calculate density-dependent mortality
            push!(deaths, o)
        else
            orgs[o].age += 1
        end
    end
    
    #unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println(sim, "Dying orgs: $(length(deaths))")
        #println("Dying dens-indep orgs: $(length(deaths))")
    end
    
    deleteat!(orgs, deaths) #delete the ones that are already dying due to mortality rate, so that they can´t also die due to density-dependent

    #Density-dependent mortality
    K = 2*(3.5/100)*(length(find(x -> x == true, landavail))*25) #7 tons/ha = 7g/100

    while currentk >= K
        println("Current biomass bigger than capacity.")

        seeds = find(x -> x.stage == "s", orgs)
        juvs = find(x -> x.stage == "j", orgs)
        adts = find(x -> x.stage == "a", orgs)
        
        if length(juvs) >= 1  # Juveniles die first         
            o = rand(juvs)
            over = /(currentk - K,K)
            Bm = orgs[o].b0am * (sum(collect(values(orgs[o].mass))))^(-1/4)*exp(-aE/(Boltz*T))
            mprob = 1 - exp(-Bm) + over

            # fix mortality probability values 
            if mprob < 0
                mprob = 0
            elseif mprob > 1
                mprob = 1
            end
            
            if 1 == rand(Distributions.Binomial(1,mprob),1)[1] || sum(values(orgs[o].mass)) <= 0
                currentk -= sum(values(orgs[o].mass)) #update currentk to calculate density-dependent mortality
                deleteat!(orgs,o)
            else
                orgs[o].age += 1
            end
        elseif length(adts) >= 1 # if there are none, adults second
            o = rand(adts)
            over = /(currentk - K,K)
            Bm = orgs[o].b0am * (sum(collect(values(orgs[o].mass))))^(-1/4)*exp(-aE/(Boltz*T))
            mprob = 1 - exp(-Bm) + over

            # fix mortality probability values 
            if mprob < 0
                mprob = 0
            elseif mprob > 1
                mprob = 1
            end
            
            if 1 == rand(Distributions.Binomial(1,mprob),1)[1] || sum(values(orgs[o].mass)) <= 0
                currentk -= sum(values(orgs[o].mass)) #update currentk to calculate density-dependent mortality
                deleteat!(orgs,o)
            else
                orgs[o].age += 1
            end
        else
            # Seeds at last
            o = rand(seeds)
            over = /(currentk - K,K)
            Bm = orgs[o].b0em * (sum(collect(values(orgs[o].mass))))^(-1/4)*exp(-aE/(Boltz*T))
            mprob = 1 - exp(-Bm) + over

            # fix mortality probability values 
            if mprob < 0
                mprob = 0
            elseif mprob > 1
                mprob = 1
            end
            
            if 1 == rand(Distributions.Binomial(1,mprob),1)[1] || sum(values(orgs[o].mass)) <= 0
                currentk -= sum(values(orgs[o].mass)) #update currentk to calculate density-dependent mortality
                deleteat!(orgs,o)
            else
                orgs[o].age += 1
            end
        end  
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
Kill organisms that where in the lost habitats.
"""
function destroyorgs!(orgs::Array{Organisms.Organism,1}, landavail::Array{Bool,2}, settings::Dict{String,Any})
    #KILL ORGANISMS in the destroyed
    kills = []
    for o in 1:length(orgs)
        x,y,f = orgs[o].location[1],orgs[o].location[2],orgs[o].location[3]
        if landavail[x,y,f] == false
            push!(kills,o)
        end        
    end
    if length(kills) > 0 # delete at index 0 generates an error
        deleteat!(orgs, kills)
    end
    #unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println("Killed orgs: $(length(kills))")
    end
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
