"""
This file contains data structures used to represent entities in the model and to store values somehow related to them
"""
# Data structures and submodel functions related to the Plant entity
# -----------------------------------------------------------------

# Data structure holding individual trait values
mutable struct Plant
    id::String
    stage::String #e,j,a
    location::Tuple # (x,y)
    sp::String #sp id, easier to read
    kernel::String
    clonality::Bool
    #### Evolvable traits ####
    compartsize::Float64
    span::Int64
    firstflower::Int64
    floron::Int64
    floroff::Int64
    seednumber::Int64
    seedon::Int64
    seedoff::Int64
    bankduration::Int64
    #### State variables ####
    age::Int64 # control death when older than max. lifespan
    mass::Dict{String, Float64}
    mated::Bool
end

# Data structures and submodel functions related to the Land entity
# -----------------------------------------------------------------

#Store simulation parameters referent to the landscape
mutable struct LandPars
    # Original landscape parameters
    npatches::Int # patches of habitat (fragments are after disturbance)
    plength::Array{Int64,1}
    initialarea::Int64
    initialconnect::Any # might not be simulated.
    # Disturbed landscape parameters
    nfrags::Any # might not be simulated
    flength::Any # might not be simulated
    disturbarea::Any # might not be simulated
    disturbconnect::Any # might not be simulated
    # general
    bufferarea::Any # might not be simulated
    meantempts::Array{Float64,1} # all fragments get the same temperature
end

mutable struct NeutralLandPars
    initialland::Array{Int64,2}
    disturbland
    meantempts::Array{Float64,1}
end