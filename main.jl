#!/usr/bin/env julia

# Get model directory and include it in Julia's loading path
cd("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/") # Only in ATOM
EDDir = pwd()
push!(LOAD_PATH,EDDir)

srand(123)

# Load Julia & model packages
using Distributions
using JLD
using Setworld
using Fileprep
using Organisms

const Boltz = 8.617e-5 # Boltzmann constant eV/K (non-SI) 1.38064852e-23 J/K if SI
const aE = 0.69 # activation energy kJ/mol (non-SI), 0.63eV (MTE - Brown et al. 2004)
const plants_gb0 = exp(25.2) # plant biomass production (Ernest et al. 2003) #TODO try a way of feeding those according to funcitonal group
const plants_fb0 = exp(26.0) # fertility rate
const tK = 273.15 # °C to K converter
const plants_mb0 = exp(19.2)

"""
    read_initials(simparams, initorgs)
Reads in and stores landscape conditions and organisms from `"landscape_init.in"` and `"organisms.in"` and stores values in composite types.
"""
function read_initials()
    #TODO dictionnary?
    simparams = Setworld.Simpars()
    begin
        env = open("landscape_init.in")
        readline(env); simparams.fxlength = tuple(parse(Int64,readline(env)))
        readline(env); simparams.fylength = tuple(parse(Int64,readline(env)))
        readline(env); simparams.fmeantemp = tuple(parse(Float64,readline(env)))
        readline(env); simparams.ftempsd = tuple(parse(Float64,readline(env)))
        readline(env); simparams.fmeanprec = tuple(parse(Float64,readline(env)))
        readline(env); simparams.fprecsd = tuple(parse(Float64,readline(env)))
        readline(env); simparams.nfrags = parse(Int64,readline(env))
        close(env)
    end
    #verify that
    simparams.nfrags == length(simparams.fxlength)

    initorgs = Organisms.InitOrgs()
    begin
        #TODO check the organisms file format: so far, all fragments get the same sps
        orgf = open("organisms.in")
        readline(orgf); initorgs.fgroups = tuple(readline(orgf)) #doest need parse for string
        readline(orgf); initorgs.sps = tuple(readline(orgf))
        readline(orgf); initorgs.init_stage = tuple(readline(orgf))
        readline(orgf); initorgs.init_abund = tuple(parse(Int64,readline(orgf)))
        #readline(orgf); genotypes = (readline(orgf))
        readline(orgf); initorgs.biomassμ = tuple(parse(Float64,readline(orgf)))
        readline(orgf); initorgs.biomasssd = tuple(parse(Float64,readline(orgf)))
        readline(orgf); initorgs.dispμ = tuple(parse.(Float64,readline(orgf)))
        readline(orgf); initorgs.dispshp = tuple(parse.(Float64,readline(orgf)))
        readline(orgf); initorgs.radius = tuple(parse(Int64,readline(orgf)))
        #check why if parsed into tuple, becomes a float
        close(orgf)
    end
    return simparams, initorgs
end

# function read_landscape()
#     #read and store initial conditions in dictionnary
#     # for nas linhas do input
#     file = open(readcsv, "landscpinit.csv")
#     map(x,y -> Dict(x => y), file(1,1:end),file(2:end,1:end))
#     landpars = Dict{Any,Float64}(file(1,1:end),file(2:end,1:end))
#
#     #Alternative
#     file = readlines("landscpinit.csv")
#     landpars = Dict()
#     simpars = Simpars()
#     for l in 1:length(file)
#         simpars.???
#         \
#         parse(str, start; greedy=true, raise=true)
#     end
# end

"""
    simulate!()
"""
function simulate()
#   INITIALIZATION
    simparams, initorgs = read_initials()

    mylandscape = landscape_init(simparams)
    orgs = newOrgs(mylandscape, initorgs)

    # INITIALIZATION ": multidimensional
    #read_landscape()
    #read_orgs()

# MODEL RUN
    for t in 1:timesteps
        #competition on 2 steps: reflects on compterm
        #develop!()
        projvegmass!(mylandscape,orgs)
        nogrowth = allocate!(mylandscape,orgs,t,aE,Boltz)
        #TODO check if there is no better way to keep track of individuals that are not growing
        # offspring = reproduce(mylandscape,orgs)
        reproduce!(mylandscape,orgs)
        if rem(timestep - 22, 52) == 0
            disperse!(orgs,mylandscape)
        end
        establish!(mylandscape, orgs)
        survive!(mylandscape, orgs,nogrowth)
        # Disturbances:
        ## Dynamical landscape change
        # if t #something
        #     function update_landscape!()
        # end
        ## Invasion
        # read_orgs(invasivefile)

        # Output:
        #orgs
        outputorgs()
        save(string("week",t))
        #network interactions
        #outputnetworks()
        #save(string("week",4*t)) more reasonable interval
    #end
    return mylandscape, orgs_init
end

simulate()

"""
    outputorgs()
Saves a long format table (`.tsv` file) with the organisms field informations.
"""
function outputorgs(orgs)

    #mk dir with simulation parameters identifier

    sep = "\t"

    output = open(string("orgsweek",t,".tsv"), "w")
    print(output, "id", sep)
    print(output, "location", sep)
    print(output, "species", sep)
    print(output, "stage", sep)
    print(output, "functional_group", sep)
    print(output, "genotype", sep)
    print(output, "dispersal_pars", sep)
    println(output)

    for o in 1:length(orgs)
        writedlm(output, [orgs[o].id orgs[o].location orgs[o].sp orgs[o].stage orgs[o].fgroup orgs[o].genotype orgs[o].disp])
    end

    close(output)

end


#TODO Store organisms parameters that regulate model run: Might interfere with the initialization values??
mutable struct orgpars
end

#Dictionnary stores functio2nal groups parameters
"""
    fgpars()
Stores functional groups parameters in a dictionnary that is consulted for every function involving organism's simulation.
    !!!! Better than struct in this case because it is possible to write general funcitons that match the strings identifying the fg og the organism and the key in the dictionnary
"""
# function fgparsdict()
#     fgpars = Dict()
#     # parse files in EDDir/functionalgroups and for each file, create an entry in the file with the first 3 letters of the group and the parametr it controls. takes the parameters of a list
#     return fgpars
# end
