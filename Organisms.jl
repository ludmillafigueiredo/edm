
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

export Organism, OrgsRef, initorgs, projvegmass!, compete, develop!, allocate!, checkboundaries, reproduce!, mate!, mkoffspring!, disperse!, germinate, establish!, survive!, shedd!, destroyorgs!, release!

# Set up model constants
#const Boltz = 1.38064852e-23 # Alternatively:
const Boltz = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
#const aE = 1e-19 # Alternatively:
const aE = 0.63 #0.65 - eV Brown & Sibly MTE book chap 2
const µ_wind = 0.1
const λ_wind = 3
const µ_ant = 1
const λ_ant = 0.2

# Initial organisms parametrization is read from an input file and stored in OrgsRef
mutable struct OrgsRef
sp_id::Array{String, 1}
kernel::Dict{String,String}
kernel_sd::Dict{String,Float64}
e_mu::Dict{String,Float64}
e_sd::Dict{String,Float64}
seedbank::Dict{String,Int}
b0g::Dict{String,Float64}
b0g_sd::Dict{String,Float64}
b0em::Dict{String,Float64}
b0em_sd::Dict{String,Float64}
b0jm::Dict{String,Float64}
b0jm_sd::Dict{String,Float64}
b0am::Dict{String,Float64}
b0am_sd::Dict{String,Float64}
b0jg::Dict{String,Float64}
b0jg_sd::Dict{String,Float64}
b0ag::Dict{String,Float64}
b0ag_sd::Dict{String,Float64}
sestra::Dict{String,Bool}
#dyad::Dict{String,Float64}
max_seedn::Dict{String,Int}
floron::Dict{String,Int}
floron_sd::Dict{String,Int}
floroff::Dict{String,Int}
floroff_sd::Dict{String,Int}
seedon::Dict{String,Int}
seedon_sd::Dict{String,Int}
seedoff::Dict{String,Int}
seedoff_sd::Dict{String,Int}
max_mass::Dict{String,Float64}
max_span::Dict{String,Int}
max_span_sd::Dict{String,Int}
mass_mu::Dict{String,Float64}
mass_sd::Dict{String,Float64}
abund::Dict{String,Int}
end

mutable struct Organism
id::String
stage::String #e,j,a
location::Tuple # (x,y,frag)
sp::String #sp id, easier to read
mass::Dict
#### Evolvable traits ####
kernel::String
e_mu::Float64
seedbank::Int64
b0g::Float64
b0em::Float64
b0jm::Float64
b0am::Float64
b0jg::Float64
b0ag::Float64
sestra::Bool
#dyad::Float64
floron::Int
floroff::Int
wseedn::Int
seedon::Int
seedoff::Int
max_mass::Float64
max_span::Int64
#topt::Float64
#### State variables #### #TODO put all state variable together
age::Int64 # control death when older than max. lifespan
mated::Bool
end
Organism(id,stage,location,sp,mass,kernel,e_mu,seedbank,b0g,b0em,b0jm,b0am,b0jg,b0ag,sestra,floron,floroff,wseedn,seedon,seedoff,max_mass,max_span) = Organism(id,stage,location,sp,mass,kernel,e_mu,seedbank,b0g,b0em,b0jm,b0am,b0jg,b0ag,sestra,floron,floroff,wseedn,seedon,seedoff,max_mass,max_span,26,false)

"""
    initorgs(landavail, orgsref,id_counter)

Initializes the organisms characterized in the input info stored in `orgsref` and distributes them in the available landscape `landavail`. Stores theindividuals in the `orgs` array, which holds all organisms being simulated at any given time.

"""
function initorgs(landavail::Array{Bool,N} where N, orgsref::Organisms.OrgsRef, id_counter::Int)

orgs = Organism[]

for frag in 1:size(landavail,3)

for s in orgsref.sp_id # all fragments are populated from the same species pool

XYs = hcat(rand(1:size(landavail,1),orgsref.abund[s]),
rand(1:size(landavail,2),orgsref.abund[s]))

for i in 1:orgsref.abund[s]

id_counter += 1 # update individual counter

neworg = Organism(hex(id_counter),
rand(["a" "j" "e"]),
(XYs[i,1],XYs[i,2],frag),
s,
Dict("veg" => 0.0,
"repr" => 0.0),
orgsref.kernel[s],
orgsref.e_mu[s],
orgsref.seedbank[s],
rand(Distributions.Normal(mean(Uniform(min(orgsref.b0g[s],orgsref.b0g_sd[s]),
max(orgsref.b0g[s],orgsref.b0g_sd[s])+0.00000001)),
abs(-(orgsref.b0g[s],orgsref.b0g_sd[s])/6))), #b0g
rand(Distributions.Normal(mean(Uniform(min(orgsref.b0em[s],orgsref.b0em_sd[s]),
max(orgsref.b0em[s],orgsref.b0em_sd[s])+0.00000001)),
abs(-(orgsref.b0em[s],orgsref.b0em_sd[s])/6))), #b0em
rand(Distributions.Normal(mean(Uniform(min(orgsref.b0jm[s],orgsref.b0jm_sd[s]),
max(orgsref.b0jm[s],orgsref.b0jm_sd[s])+0.00000001)),
abs(-(orgsref.b0jm[s],orgsref.b0jm_sd[s])/6))), #b0jm
rand(Distributions.Normal(mean(Uniform(min(orgsref.b0am[s],orgsref.b0am_sd[s]),
max(orgsref.b0am[s],orgsref.b0am_sd[s])+0.00000001)),
abs(-(orgsref.b0am[s],orgsref.b0am_sd[s])/6))), #b0am
rand(Distributions.Normal(mean(Uniform(min(orgsref.b0jg[s],orgsref.b0jg_sd[s]),
max(orgsref.b0jg[s],orgsref.b0jg_sd[s])+0.00000001)),
abs(-(orgsref.b0jg[s],orgsref.b0jg_sd[s])/6))), #b0jg
rand(Distributions.Normal(mean(Uniform(min(orgsref.b0ag[s],orgsref.b0ag_sd[s]),
max(orgsref.b0ag[s],orgsref.b0ag_sd[s])+0.00000001)),
abs(-(orgsref.b0ag[s],orgsref.b0ag_sd[s])/6))), #b0ag
orgsref.sestra[s],
Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.floron[s],orgsref.floron_sd[s]),
max(orgsref.floron[s],orgsref.floron_sd[s])+0.00000001)),
abs(-(orgsref.floron[s], orgsref.floron_sd[s])/6))),RoundUp)), #floron
Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.floroff[s], orgsref.floroff_sd[s]),
max(orgsref.floroff[s], orgsref.floroff_sd[s])+0.00000001)),
abs(-(orgsref.floroff[s], orgsref.floroff_sd[s])/6))),RoundUp)), #floroff
0,
Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.seedon[s], orgsref.seedon_sd[s]),
max(orgsref.seedon[s], orgsref.seedon_sd[s])+0.00000001)),
abs(-(orgsref.seedon[s],orgsref.seedon_sd[s])/6))),RoundUp)), #seedon
Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.seedoff[s], orgsref.seedoff_sd[s]),
max(orgsref.seedoff[s], orgsref.seedoff_sd[s])+0.00000001)),
abs(-(orgsref.seedoff[s],orgsref.seedoff_sd[s])/6))),RoundUp)), #seedoff
0.0, #max_mass
Int(round(rand(Distributions.Normal(mean(Uniform(min(orgsref.max_span[s], orgsref.max_span_sd[s]),
max(orgsref.max_span[s], orgsref.max_span_sd[s])+0.00000001)),
abs(-(orgsref.max_span[s], orgsref.max_span_sd[s])/6))),RoundUp)))#max_span

# adult max biomass
if neworg.e_mu == 0.0001
neworg.max_mass = 1
elseif neworg.e_mu == 0.0003
neworg.max_mass = 2
elseif neworg.e_mu == 0.001
neworg.max_mass = 5
else
error("Check seed sizes in input")
end

# weekly number of seeds
neworg.wseedn = Int(round(orgsref.max_seedn[s]/(neworg.floroff-neworg.floron+1)))

# initial biomass
if neworg.stage in ["e","j"]
neworg.mass["veg"] = neworg.e_mu
else
neworg.mass["veg"] = neworg.max_mass * 0.5
end
push!(orgs, neworg)

end
end
end

return orgs, id_counter

end

"""
allocate!(landscape, orgs, t, aE, Boltz, setting, orgsref, T)
Calculates biomass gain according to the metabolic theory (`aE`, `Boltz` and `T` are necessary then). According to the week being simulated, `t` and the current state of the individual growing ( the biomass gained is
"""
function allocate!(landscape::Array{Dict{String, Float64},N} where N, orgs::Array{Organism,1}, t::Int64, aE::Float64, Boltz::Float64, settings::Dict{String, Any},orgsref::Organisms.OrgsRef,T::Float64)
#1. Initialize storage of those that dont growi and will have higher prob of dying (later)
nogrowth = Int64[]

growing = find(x->(x.stage == "a" || x.stage == "j"),orgs)
for o in growing

# we verified that a species lifespan determines how fast it should grow
if orgs[o].max_span <= 1
b0 = orgs[o].b0ag*100
elseif 1 < orgs[o].max_span <= 2
b0 = orgs[o].b0ag*15
elseif orgs[o].max_span > 2
b0 = orgs[o].b0ag*10
else
error("Check species maximum lifespan.")
end

#only vegetative biomass helps growth
grown_mass = b0*(orgs[o].mass["veg"])^(3/4)*exp(-aE/(Boltz*T))

if isapprox(grown_mass,0) # if it is not growing, there is no need to allocate
push!(nogrowth,o)
elseif orgs[o].stage == "j"
# juveniles grow vegetative biomass only
orgs[o].mass["veg"] += grown_mass 
elseif orgs[o].stage == "a" &&
(orgs[o].floron <= rem(t,52) < orgs[o].floroff) &&
(sum(collect(values(orgs[o].mass))) >= 0.5*(orgs[o].max_mass))
# adults in their reproductive season and with enough weight, invest in reproduction
#sowingmass = (5.5*(10.0^(-2)))*((orgs[o].mass["veg"]/(orgs[o].floroff-orgs[o].floron + 1) + grown_mass)^0.95)
if haskey(orgs[o].mass,"repr")
orgs[o].mass["repr"] += grown_mass #sowingmass 
else
orgs[o].mass["repr"] = grown_mass #sowingmass
end
elseif orgs[o].stage == "a" && orgs[o].mass["veg"] < orgs[o].max_mass
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
if  orgs[j].mass["veg"] >= 0.5*orgs[j].max_mass
# If an individual grows quite fast, it is more vigorous, and should transfer it to adult fecundity. The only variable capable of transfering this property is the weigh, which, combined with the MTE rate, makes it  generate more offspring
orgs[j].stage = "a"
end
end
end

"""
checkboundaries(sourcefrag,xdest, ydest, fdest)
`source` and `dest` contain the location indexes of the source (mother plant) and the pollen/seed. `checkboundaires()` verifies whether the new polen/seed location `(x,y)` is inside a habitat fragment (same as the source -`frag`- or another one insed the patch). Return a boolean that controls whether the process (reproduction or emergency/germination) proceeds or not.
"""
function checkboundaries(landavail::Array{Bool,N} where N, xdest::Int64, ydest::Int64, fdest::Int64) #TODO is this function necessary?
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
mate!()
Calculate proportion of insects that reproduced (encounter?) and mark that proportion of the population with the `mated` label.
- visited: reduction in pollination service
"""
function mate!(orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any}, scen::String, regime:: String, td::Int, visited::Int)

ready = find(x-> x.stage == "a" && x.mass["repr"] > x.e_mu, orgs) # TODO find those with higher reproductive mas than the mean nb of seeds * seed mass.
pollinated = []
npoll = 0

if length(ready) > 0 # check if there is anyone flowering
# Scenarios were pollination is not species-specific
if scen in ["indep" "equal"]  # calculation of number of pollinated individuals is different, but the actual pollination (non species-specific) is

# Determine number of individuals that get pollinated (species is not relevant)
if scen == "indep"
npoll = Int(ceil(rand(Distributions.Uniform(0.5, 1))[1] * length(ready) * 1))
println("Scenario of INDEP pollination loss")

elseif scen == "equal" #all species lose pollination randomly (not species-specific)
println("Scenario of EQUAL pollination loss")
# Check the regime
if t >= td 
if regime == "pulse"
if t != td # t is not td more often, therefore it is simpler to avoid checking for the less simpler situation most of the time 
npoll = Int(ceil(rand(Distributions.Uniform(0.5,1))[1] * length(ready) * 1))
else
npoll = Int(ceil(rand(Distributions.Uniform(0.5,1))[1] * length(ready) * visited))
end

elseif regime == "exp"

npoll = Int(ceil(rand(Distributions.Uniform(0.5,1))[1] * length(ready) * exp(-0.5*(t-td))))
# exp(tp - t) makes the pollination loss decrease from 1 (tp = t) to 0
println("Regime of EXP pollination loss, with $(length(ready)) being ready and $npoll being pollinated")

elseif regime == "const"
#if t <= tdend 
#    npoll = Int(ceil(rand(Distributions.Uniform(0.5,1))[1] * length(ready) * visited))
#else
npoll = Int(ceil(rand(Distributions.Uniform(0.5,1))[1] * length(ready) * 1)) # pollination goes back to 'normal' after press disturbance
#end
else
error("Please choose a pollination regime \"regime\" in insect.csv:
- \"pulse\":
- \"const\":
- \"exp\"")                    
end
else
npoll = Int(ceil(rand(Distributions.Uniform(0.5,1))[1] * length(ready) * 1))
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
function mkoffspring!(orgs::Array{Organisms.Organism,1}, t::Int64, settings::Dict{String, Any},orgsref::Organisms.OrgsRef, id_counter::Int)

offspring = Organism[]

# Asexually produced offspring
asexuals = filter(x -> x.mated == false && x.sestra == true && x.mass["repr"] > x.e_mu, orgs) #find which the individuals the can reproduce assexually and then go through them, by species

for sp in unique(getfield.(asexuals, :sp))
cloning = find(x -> x.sp == sp && x.id in unique(getfield.(asexuals, :id)) , orgs)

# start production counting
spclonescounter = 0

for c in cloning # mothers cloning
offs = div(0.5*orgs[c].mass["repr"], orgs[c].e_mu)

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
orgs[c].mass["repr"] -= (offs * orgs[c].e_mu)

# get a copy of the mother, which the clones will look like 
clonetemplate = deepcopy(orgs[c])
clonetemplate.stage = "j" #clones have already germinated
clonetemplate.mass["veg"] = orgs[c].e_mu
clonetemplate.mass["repr"] = 0.0

for o in offs

id_counter += 1

clone = deepcopy(clonetemplate)
# clone's own traits
clone.id = hex(id_counter)
clone.location = (clonetemplate.location[1] + rand(Distributions.Bernoulli())[1],
clonetemplate.location[2] + rand(Distributions.Bernoulli())[1],
clonetemplate.location[3]) # clones are spread in one of the neighboring cells - or in the same as the mother
push!(offspring, clone)
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

emu = orgs[s].e_mu
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

#newvalue = rand(Distributions.Normal(0,abs(embryo.e_mu-conspp.e_mu)/embryo.e_mu))
#embryo.e_mu + newvalue >= orgsref.e_mu[embryo.sp] ? # if seed biomass or minimal biomass would smaller than zero, it does not chenge
#embryo.e_mu += newvalue : embryo.e_mu += 0
embryo.b0g += rand(Distributions.Normal(0,abs(embryo.b0g-conspp.b0g+0.0000001)/embryo.b0g))
embryo.b0em += rand(Distributions.Normal(0,abs(embryo.b0em-conspp.b0em+0.0000001)/embryo.b0em))
embryo.b0jm += rand(Distributions.Normal(0,abs(embryo.b0jm-conspp.b0jm+0.0000001)/embryo.b0jm))
embryo.b0am += rand(Distributions.Normal(0,abs(embryo.b0am-conspp.b0am+0.0000001)/embryo.b0am))
embryo.b0jg += rand(Distributions.Normal(0,abs(embryo.b0jg-conspp.b0jg+0.0000001)/embryo.b0jg))
embryo.b0ag += rand(Distributions.Normal(0,abs(embryo.b0ag-conspp.b0ag+0.0000001)/embryo.b0ag))
embryo.floron += Int(round(rand(Distributions.Normal(0,abs(embryo.floron-conspp.floron+0.0000001)/embryo.floron)),RoundUp))
embryo.floroff += Int(round(rand(Distributions.Normal(0,abs(embryo.floroff-conspp.floroff+0.0000001)/conspp.floroff)),RoundUp))
embryo.seedon += Int(round(rand(Distributions.Normal(0,abs(embryo.seedon-conspp.seedon+0.0000001)/embryo.seedon)),RoundUp))
embryo.seedoff += Int(round(rand(Distributions.Normal(0,abs(embryo.seedoff-conspp.seedoff+0.0000001)/embryo.seedoff)),RoundUp))
#newvalue = Int(round(rand(Distributions.Normal(0,abs(embryo.max_span-conspp.max_span/embryo.max_span)),RoundUp))
#embryo.max_span + newvalue > orgsref.max_span[embryo.sp] ?
#embryo.max_span += newvalue : embryo.max_span += 0

# set embryos own individual non-evolutionary traits
embryo.id = hex(id_counter) 
embryo.mass = Dict("veg" => embryo.e_mu,
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
seedsi = find(x -> x.stage == "e" && x.age == 0 && x.seedon <= rem(t,52) < x.seedoff, orgs) # using a condition "outside" orgs might not work. This condition with orgsref only works because orgsref always has the sp names of x as keys in the dictionnary. If presented with a key that it does do contain, it throws an error.
return seedsi
end

"""
disperse!(landscape, orgs, t, seetings, orgsref,)
Seeds are dispersed.
"""
function disperse!(landavail::Array{Bool,N} where N,seedsi, orgs::Array{Organisms.Organism, 1}, t::Int, settings::Dict{String, Any}, orgsref::Organisms.OrgsRef, landpars::Any, tdist::Any)#Setworld.LandPars)

lost = Int64[]

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

# Find patch
θ = rand([0 0.5π π 1.5π]) # TODO tentar rand(0:2)*pi?
xdest = orgs[d].location[1] + dist*round(Int64, cos(θ), RoundNearestTiesAway)
ydest = orgs[d].location[2] + dist*round(Int64, sin(θ), RoundNearestTiesAway)
fsource = orgs[d].location[3] 

if checkboundaries(landavail, xdest, ydest, fsource) && landavail[xdest, ydest, fsource] == true # intra fragment dispersal
orgs[d].location = (xdest,ydest,fsource)

else # check if there are other fragments that can be reached
if landpars.initialconnect != nothing || landpars.disturbconnect != nothing
# check if    
if settings["disturbtype"] in ["none" "temp" "poll"]
Ah = landpars.initialarea
connects = landpars.initialconnect
else
if t in tdist
Ah = landpars.disturbarea
connects = landpars.disturbconnect
else
Ah = landpars.initialarea
connects = landpars.initialconnect
end
end

AT = landpars.bufferarea

# check if the fragment is reached
fdest = find(x -> (x <= dist && x > 0), connects[:,fsource]) |> (y -> (length(y) > 0 ? fdest = rand(x) : false))
else
fdest = false
end

#reset location or kill individual
if fdest == false
push!(lost,d)
else
rand(Distributions.Bernoullli(1,Ah/AT))[1] == 1 ? orgs[d].location = (xdest,ydest,fdest) : push!(lost,d)
end
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
Bg = org.b0g * (org.mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
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
Seeds only have a chance of establishing in patches not already occupied by the same funcitonal group, in. When they land in such place, they have a chance of germinating (become seedlings - `j` - simulated by `germinate!`). Seeds that don't germinate stay in the seedbank, while the ones that are older than one year are eliminated.
"""
function establish!(landscape::Array{Dict{String,Float64},N} where N, orgs::Array{Organisms.Organism,1}, t::Int, settings::Dict{String, Any}, orgsref::Organisms.OrgsRef, T::Float64)
#REFERENCE: May et al. 2009
establishing = find(x -> x.stage == "e", orgs)

lost = Int64[]

for o in establishing
if orgs[o].seedon <= rem(t,52) #only after release and dispersal a seed can establish
orgcell = orgs[o].location
# TODO: take it out, it is obolete, since there is no more biomass projection 

if germinate(orgs[o],T,settings)
orgs[o].stage = "j"
#unity test
#open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
#    println(sim, "Became juvenile: $(orgs[o].id)")
#end
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
function survive!(orgs::Array{Organisms.Organism,1}, t::Int, cK::Float64, K::Float64, settings::Dict{String, Any}, orgsref::Organisms.OrgsRef, landavail::Array{Bool,N} where N,T, nogrowth::Array{Int64,1})
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
if orgs[o].age >= orgs[o].seedbank
mprob = 1
elseif rem(t,52) > orgs[o].seedoff #seeds that are still in the mother plant cant die. If their release season is over, it is certain thatthey are not anymore, even if they have not germinated 
Bm = orgs[o].b0em * (orgs[o].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
#println("Bm: $Bm, b0em = $(orgs[o].b0em), seed mass = $(orgs[o].mass["veg"])")
mprob = 1 - exp(-Bm)
else
mprob = 0
end 
elseif (orgs[o].age/52) >= orgs[o].max_span*52 #oldies die
mprob = 1
elseif orgs[o].stage == "j"
Bm = orgs[o].b0jm * (orgs[o].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
mprob = 1 - exp(-Bm)
else #adults
Bm = orgs[o].b0am * (orgs[o].mass["veg"]^(-1/4))*exp(-aE/(Boltz*T))
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

#unity test
#open(string("EDoutputs/",settings["simID"],"/simulog.txt"),"a") do sim
#    println(sim,"Dying dens-indep orgs: $(length(deaths))")
#end
deleteat!(orgs, deaths) #delete the ones that are already dying due to mortality rate, so that they can´t also die due to density-dependent
open(abspath(joinpath(settings["outputat"],settings["simID"],"simulog.txt")),"a") do sim
println(sim, "Density-independent mortality: ", length(deaths))
end
deaths = Int64[] # reset

## Density-dependent mortality
if sum(map(x -> x.mass["veg"],orgs)) > K
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
if length(fullcells) > 0
for c in fullcells
#find those that are in the same grid
samegrid = filter(x -> x.location == locs[c], orgs)
# unity test
#open(string("EDoutputs/",settings["simID"],"/simulog.txt"), "a") do sim
#    println(sim,"N of inds in same cell: $(length(samegrid))")
#end
#println("Inds in same cell $(locs[c]) : $samegrid")
# sum their weight to see if > than carrying capacity.
while sum(map(x -> x.mass["veg"],samegrid)) > cK
# any individual can die
d = rand(samegrid,1)[1]

if d.stage == "a"
Bm = d.b0am * (sum(collect(values(d.mass["veg"]))))^(-1/4)*exp(-aE/(Boltz*T))
mprob = 1 - exp(-Bm)
elseif d.stage == "j"
Bm = d.b0jm * (sum(collect(values(d.mass["veg"]))))^(-1/4)*exp(-aE/(Boltz*T))
mprob = 1 - exp(-Bm)
else
Bm = d.b0em * (sum(collect(values(d.mass["veg"]))))^(-1/4)*exp(-aE/(Boltz*T))
mprob = 1 - exp(-Bm)
end

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
samegrid = filter(x -> x.location == locs[i], orgs)
#open(string("EDoutputs/",settings["simID"],"simulog.txt"), "a") do sim
#    println("$(orgs[o].stage) dying DEP.")
#end
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
Kill organisms that where in the lost habitats.
"""
function destroyorgs!(orgs::Array{Organisms.Organism,1}, landavail::Array{Bool,N} where N, settings::Dict{String,Any})
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
