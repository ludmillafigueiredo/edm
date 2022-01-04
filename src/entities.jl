"""
This file contains data structures used to represent entities in the model and to store values somehow related to them
"""

"""
using DataFrames
"""

# Data structures and submodel functions related to the Plant entity
# -----------------------------------------------------------------

# Data structure holding plant mass values
mutable struct Mass
    repr::Float64
    root::Float64
    stem::Float64
    leaves::Float64
end

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
    mass::Mass
    mated::Bool
    flagged::Bool
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

#DataStructure holding the population of a Cell in the landscape (with some properties...maybe)
# struct Cell
#     plants::Vector{Plant}
#     habitability::Bool
#     location::Tuple
# end

struct Landscape
    plants::Array{Vector{Plant}, 2}
    celllocks::Array{ReentrantLock, 2}
    dispersal::Array{Vector{Plant}, 2}
    habitability::BitArray{2}
end

#DataStructure holding the settings
struct Settings
    output_freq::Int64
    temp_ts::String
    pollination::String
    disturb_type::String
    outputat::String
    initial_land::String
    sppinput::String
    disturb_land::Any
    rseed::Int64
    simID::String
end
