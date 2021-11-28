"""
This file contains data structures used to represent entities in the model and to store values somehow related to them
"""

"""
using DataFrames
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
    pollen_vector::String
    self_failoutcross::Bool
    self_proba::Float64
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

# Data structure for storing pollination
mutable struct PollPars
    scen::String
    regime::Any # it can be Nothing or a DataFrames.DataFrame
end


# Data structures related to the Land entity
# ------------------------------------------

#Store files describing the landscape
mutable struct LandPars
    initial::String
    disturbance::Any # it can be Nothing or a DataFrames.DataFrame
end
