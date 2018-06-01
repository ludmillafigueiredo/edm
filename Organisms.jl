"""
This module contains the

Organisms have the same attributes, whose specific values differ according to functional groups (or not?). They interact when in the vicinity of each other (this might be detected over a certain distance or not - change the range of search).
"""
module Organisms

using Distributions
using JuliaDB
using DataValues
using Setworld
using Fileprep

export Organism, OrgsRef, newOrgs, projvegmass!, compete, develop!, allocate!, meanExP, checkboundaries, reproduce!, disperse!, germinate, establish!, survive!

#TODO put them in OrgsRef
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
kernel::Dict{String,Int}
e_mu::Dict{String,Float64}
e_sd::Dict{String,Float64}
b0g::Dict{String,Float64}
b0em::Dict{String,Float64}
b0am::Dict{String,Float64}
b0eg::Dict{String,Float64}
b0ag::Dict{String,Float64}
sestra::Dict{String,Int}
dyad::Dict{String,Float64}
floron::Dict{String,Int}
floroff::Dict{String,Int}
sripe::Dict{String,Int}
seedon::Dict{String,Int}
seedoff::Dict{String,Int}
max_mass::Dict{String,Float64}
max_span::Dict{String,Int64}
mass_mu::Dict{String,Float64}
mass_sd::Dict{String,Float64}
abund::Dict{String,Int}
end

mutable struct Organism
id::String
location::Tuple # (x,y,frag)
sp::String #sp id, easier to read
mass::Dict
#fgroup::String # Plant, insect or more precise functional group
stage::String #e,j,a
age::Int64 # controls passing stages, phenology and
mated::Bool # TODO is it necessary?
genotype::Array{String,2} #initialize separately
radius::Int64 # TODO reproductive and vegetative area of influence. Not Tuple because not
#Organism() = new()
end
Organism(id,location,sp,biomass) = Organism(id,location,sp,biomass,fgroup,"a", 0,false,["A" "A"],0) #this is individuals are initialized in the beginning of the simulationy
# TODO fgroups are sps for now

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
function newOrgs(landscape::Array{Setworld.WorldCell,N} where N,orgsref::OrgsRef)

    orgs = Organism[]

    for frag in 1:size(landscape,3)

        for s in orgsref.sp_id

            XYs = hcat(rand(1:size(landscape,1),orgsref.abund[s]),
                       rand(1:size(landscape,2),orgsref.abund[s]))

            for i in 1:orgsref.abund[s]
                neworg = Organism(string(s, "-", length(orgs) + 1), #sp_id isnt ind id!
                                  (XYs[i,1],XYs[i,2],frag),
                                  s,
                                  Dict("veg" => rand(Distributions.Normal(orgsref.mass_mu[s],orgsref.mass_sd[s])), "repr" => 0))
                push!(orgs, neworg)
            end
        end
    end
    #TODO different fgroups have different metabolic rates. Get it form OrgsRef
    #TODO or make sure that it has varied in the global scope as well rather export it?
    return orgs
end

"""
        projvegmass!(landscape,orgs)
        Rewrites the projected mass of each organisms stored in `orgs` into the `neighs` field of `landscape`. This projection means that the total biomass is divided into the square area delimited by the organism's `radius`.
        """
function projvegmass!(landscape::Array{Setworld.WorldCell, N} where N, orgs::Array{Organism,1}, settings::Dict{String,Any})
    
    competing = find(x->(x.stage == "a" || x.stage == "j"),orgs) #juveniles com ashard as adults, but have higher growth rate and lower mortality

    for o in competing
        x, y, frag = orgs[o].location
        orgs[o].radius = round(Int64, (sqrt(orgs[o].biomass["veg"]^(2/3)) - 1)/2, RoundUp)

        r = orgs[o].radius # separated for debugging

        #sp = orgs[o].sp

        projmass = /(orgs[o].mass["veg"], ((2*r+1)^2))

        # unity test
        # open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        #     println(sim, orgs[o].id," has ",orgs[o].biomass["veg"], " radius $r and projvegmass:", projmass) #ugly format to avoid risking some anoying errors that have been happening
        # end

        for j in (y-r):(y+r), i in (x-r):(x+r) #TODO usar a funcao da FON Q trabalha com quadrantes? dar mais peso para steming point?
            if !checkbounds(Bool,landscape[:,:,frag],j,i) # check boundaries: absorbing borders: the biomass is not re-divided to the amount of cells inside the fragment. What is projected outside the fragmetn is actually lost: Edge effect
                continue
            else
                if haskey(landscape[i,j,frag].neighs,"p")
                    landscape[i,j,frag].neighs["p"] += projmass
                else
                    landscape[i,j,frag].neighs["p"] = projmass
                end
            end
        end
    end
end

"""
        compete(landscape, orgs)
        Each plant in `orgs` will check neighboring cells (inside it's zone of influence radius `r`) and detected overlaying ZOIs. When positive, the proportion of 'free' plant biomass is calculated: (focus plant biomass - sum(non-focus vegetative biomass))/(focus plant biomass), and normalized to the total area of projected biomass .
        """
function compete(landscape::Array{Setworld.WorldCell, N} where N, org::Organism, settings::Dict{String,Any})

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
            # open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #     println(sim, "$(org.id) out of bound projection at $(org.location)")
            #end
        elseif haskey(landscape[i,j,frag].neighs,"p")
            #landscape[i,j,frag].neighs[fg] > 0 # check the neighborhood of same fgroup for competition
            nbsum += landscape[i,j,frag].neighs["p"] - /(org.mass["veg"],(2*r+1)^2) #sum vegetative biomass of neighbors only (exclude focus plant own biomass)
        end
    end

    compterm = /(org.mass["veg"] - nbsum, org.mass["veg"])
    # # unity test
    # open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
    #     println(sim, "$(org.id) radius = $r and compterm = $compterm")
    # end

    return compterm
end

"""
        allocate!(orgs, landscape, aE, Boltz, OrgsRef)
        Calculates biomass gain according to MTE rate and depending on competition. Competition is measured via a biomass-based index `compterm` (`compete` function). This term gives the proportion of actual biomass gain an individual has. If competition is too strong (`compterm` < 0), the individual has a higher probability of dying. If not, the biomass is allocated to growth or reproduction, according to the developmental `stage` of the organism and the season (`t`) (plants start allocating to week 12).
"""
function allocate!(landscape::Array{Setworld.WorldCell,N} where N, orgs::Array{Organism,1}, t::Int64, aE::Float64, Boltz::Float64, settings::Dict{String, Any})
    #1. Initialize storage of those that are ont growing and have higher prob of dying (later)
    nogrowth = Int64[]

    #TODO filter insects out?
    for o in 1:length(orgs)

	if orgs[o].stage == "e" #embryos dont allocate because they dont grow
	    continue
	else

            # 2.a Get local temperature
            T = landscape[orgs[o].location[1], orgs[o].location[2], orgs[o].location[3]].temp

            # This MTE rate comes from dry weights: fat storage and whatever reproductive structures too, but not maintenance explicitly
            # Any cost related to insufficient minimal biomass goes into the survival probability function
            # 2.b Check for competition

            compterm = compete(landscape, orgs[o], settings)
            # unity test
            #println(simulog, org.id," weights",org.biomass["veg"]," had $nbsum g overlap")
            open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                println(sim, orgs[o].id, "  compterm $compterm")
            end

            if compterm <= 0
                #2.c Those not growing will have higher chance of dying
                push!(nogrowth,o)
            else
                sp = orgs[o].sp 
                # seedling and adults grow at different rates
                if orgs[o].stage == "j"
                    b0 = orgsref.b0eg[sp]
                elseif orgs[o].stage == "a"
                    b0 = orgsref.b0ag[sp]
                else
                    error("Seed or egg trying to compete.")
                end
                
                grown_mass = b0*(compterm *sum(collect(values(orgs[o].biomass))))^(3/4)*exp(-aE/(Boltz*T))
                # unity test
                open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                    println(sim, "$(orgs[o].id) should grow $grown_mass")
                end

                #Resource allocation schedule
                #TODO make it more ellaborate and includde trade-offs
                if orgs[o].stage == "j"
                    # juveniles grow
                    orgs[o].biomass["veg"] += grown_mass
                    # unity test
                    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                        println(sim, "$(orgs[o].id)-$(orgs[o].stage) grew $grown_mass")
                    end
                elseif orgs[o].stage == "a" && orgsref.floron[sp] <= t < orgsref.floroff[sp]                    #TODO extend it to summer?
                    # adults invest in reproduction
                    #unity test
                    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                        println(sim, "$(orgs[o].id)-$(orgs[o].stage) is FLOWERING")
                    end
                    if haskey(orgs[o].mass,"reprd")
                        orgs[o].mass["reprd"] += grown_mass
                    else
                        orgs[o].mass["reprd"] = grown_mass
                    end
                elseif orgs[o].stage == "a" #&& orgs[o].biomass["veg"] < 50 # TODO refer it to a 50% of the species biomass #individuals that are too small dont reproduce
                    orgs[o].biomass["veg"] += grown_mass
                    # unity test
                    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                        println(sim, "$(orgs[o].id)-$(orgs[o].stage) grew VEG $grown_mass")
                    end
                else
                    error("Seed or egg trying to allocate.")
                end
            end

            #unity test
            # open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #     println(sim, "current biomass: $(orgs[o].biomass)")
            # end

            # unity test
            open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                println(sim, "Not growing $nogrowth")
            end
        end
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
checkboundaries(sourcefrag,xdest, ydest, fdest)
`source` and `dest` contain the location indexes of the source (mother plant) and the pollen/seed. `checkboundaires()` verifies whether the new polen/seed location `(x,y)` is inside a habitat fragment (same as the source -`frag`- or another one insed the patch). Return a boolean that controls whether the process (reproduction or emergency/germination) proceeds or not.
"""
function checkboundaries(landscape::Array{Setworld.WorldCell,N} where N, xdest::Int64, ydest::Int64, fdest::Int64) #TODO is this function necessary?
    #check inside frag
    if checkbounds(Bool, landscape[:,:,fdest], xdest, ydest)
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
    """
"""
    mate!()
    Calculate proportion of insects that reproduced (encounter?) and mark that proportion of the population with the `mated` label.
    GEnetic comes here
    """

"""
    mkoffspring!()
    After mating happened (marked in `reped`), calculate the amount of offspring
    """
function mkoffspring!(landscape::Array{Setworld.WorldCell, N} where N, orgs::Array{Organisms.Organism,N} where N, t::Int64, settings::Dict{String, Any},orgsref::OrgsRef)
    #TODO sort out reproduction mode (pollination or not) according to functional group
    # if # pollination depending plants
    #     pollination()
    # elseif "insect"
    #     parents_genes = mate!()
    # end

    ferts = find(x -> x.mated == true, orgs)
    
    #unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println(sim, "Reproducing: $reproducing week $t")
    end

    offspring = Organism[]

    for o in ferts

        emu = orgsref.e_mu[orgs[o].sp]
        offs =  round(Int64, /(orgs[o].biomass["reprd"],emu), RoundDown)
        orgs[o].biomass["reprd"] -= (offs * emu)

        #unity test
        open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            println(sim, "Number of offspring of ", orgs[o].id, ": ",offs)
        end

        for n in 1:offs
            #TODO check for a quicker way of creating several objects of composite-type

            embryo = Organism(string(orgs[o].sp,"-", (length(orgs) + length(offspring) + 1)) ,
                              orgs[o].location, #stays with mom until release
                              orgs[o].sp,
                              Dict("veg" => rand(Distributions.Normal(orgsref.e_mu[orgs[o].sp],orgsref.e_sd[orgs[o].sp]),1)[1]), #use seed size for the fgroup
                              "e",
                              0,
            false,
            ["A" "A"], #come from function
            0)

            push!(offspring, embryo)
	    #unitest
            open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                println(sim, "Pushed new org ", embryo, " into offpring")

            end
        end
    end
    append!(orgs, offspring)
    #unitest
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
	println(sim, "Total offspring this summer: " ,offspring)
    end
    
end

"""
    release!()
    Probably need to call disperse here, because not all "e"s in orgs are released at the same time.
    """

"""
    disperse!(offspring)
    Seeds are dispersed.
    """
function disperse!(landscape::Array{Setworld.WorldCell,N} where N,orgs::Array{Organisms.Organism, N} where N, t::Int, settings::Dict{String, Any})

    #TODO include insects
    # if achar plant (fgroup?)
    #  dipersar xovens
    # else
    #  dispersa adultos
    # end
    # find indexes in orgs if embryos that are in the right to be dispersed and which releasing season is happening at t
    dispersing = filter(x -> orgsref.seedon[x.sp] <= t < orgsref.seedoff[x.sp], filter(x -> (x.stage == "e" && x.age == 0), orgs)) #didnt mange to pipeline those, weird metho error
    #unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println(sim, "Dispersing: $dispersing")
    end

    lost = Int64[]

    for d in dispersing
        if orgsref.kernel[orgs[d].sp] == "ant"
            dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ_ant,λ_ant),1)[1])
        elseif orgsref.kernel[orgs[d].sp] == "wind"
            dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ_wind,λ_wind),1)[1])
        elseif orgsref.kernel[orgs[d].sp] == "windant"
            µ,λ = rand([[µ_wind λ_wind],[µ_ant λ_ant]])
            dist = Fileprep.lengthtocell(rand(Distributions.InverseGaussian(µ,λ),1)[1])
        else
            error("Check dispersal kernel input.")
        end

        #unity test
        open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            println(sim, orgs[d].id,"Calculated dispersal distance: $dist")
        end

        # Find patch
        θ = rand([0 0.5π π 1.5π]) # TODO tentar rand(0:2)*pi?
        xdest = orgs[d].location[1] + dist*round(Int64, cos(θ), RoundNearestTiesAway)
        ydest = orgs[d].location[2] + dist*round(Int64, sin(θ), RoundNearestTiesAway)
        fdest = orgs[d].location[3] #TODO is landing inside the same fragment as the source, for now

        if checkboundaries(landscape, xdest, ydest, fdest)
            orgs[d].location = (xdest,ydest,fdest)
            #unity test
            open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                println(sim, orgs[d].id," dispersed $dist")
            end
        else
            push!(lost,d)
            #unity test
            open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
                println(sim, orgs[d].id," dispersed $dist but died")
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
    gprob = 1 - exp(-orgsref.b0eg[org.sp])
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
function establish!(landscape::Array{Setworld.WorldCell,N} where N, orgs::Array{Organisms.Organism, N} where N, t::Int64, settings::Dict{String, Any})
    #REFERENCE: May et al. 2009
    establishing = find(x -> x.stage == "e", orgs)

    #unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println(sim, "Establishing seeds: $establishing")
    end

    lost = Int64[]

    for o in establishing
        orgcell = orgs[o].location
        orgfg= orgs[o].fgroup
        if haskey(landscape[orgcell[1], orgcell[2], orgcell[3]].neighs,orgfg)
            #push!(lost,o)
	    #unity test
   	    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        	println(sim, "$o falling into occupied cell")
    	    end
	end

        if germinate(orgs[o])
            orgs[o].stage = "j"
	    #unity test
   	    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        	println(sim, "Became juvenile: $o")
    	    end
            # the ones that dont germinate but are older than 1 year die
        end
        
    end
    deleteat!(orgs,lost)
end

"""
    survive!(orgs, nogrowth,landscape)
    Organism survival depends on total biomass, according to MTE rate. However, the proportionality constants (b_0) used depend on the cause of mortality: competition-related, where
    plants in nogrwth are subjected to two probability rates
    """
function survive!(landscape::Array{Setworld.WorldCell,N} where N,orgs::Array{Organisms.Organism,N} where N, nogrowth::Array{Int64,N} where N, settings::Dict{String, Any}, orgsref::OrgsRef)

    T = landscape[orgs[o].location[1], orgs[o].location[2], orgs[o].location[3]].temp

    deaths = Int64[]
    seeds = find(x -> orgs[o].stage == "e", orgs)    

    for o in 1:length(orgs)
	if o in seeds
            Bm = orgsref.b0em[orgs[o].sp] * (sum(collect(values(orgs[o].mass))))^(-1/4)*exp(-aE/(Boltz*T))
        elseif rem(orgs[o].age,52) >= orgsref.max_span(orgs[o].sp]
            Bm = 1
        elseif o in nogrowth
            Bm = orgsref.b0am[orgs[o].sp] * (sum(collect(values(orgs[o]mass))))^(-1/4)*exp(-aE/(Boltz*T))
            # unity test
            # open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
            #     println(sim,"$(orgs[o].id) $(orgs[o].stage) mortality rate $Bm")
            # end)
            
        end
        
        mprob = 1 - exp(-Bm)
        if 1 == rand(Distributions.Binomial(1,mprob),1)[1]
            push!(deaths, o)
        else
            orgs[o].age += 1
        end
    end
    #unity test
    open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
        println(sim, "Dying orgs: $deaths")
    end
    deleteat!(orgs, deaths)
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
