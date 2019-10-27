"""
This module contains data structures used to represent entities in the model and to store values somehow related to them
"""

module entities

export SppRef, TraitRanges, Plant, LandPars, NeutralLandPars

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

# Initial trait values is read from an input file and stored for reference in `sppref::SppRef`.
mutable struct SppRef
    sp_id::Array{String, 1}
    clonality::Dict{String,Bool}
    kernel::Dict{String,String}
    seedmass::Dict{String,Float64}
    compartsize::Dict{String,Float64}
    span_min::Dict{String,Float64}
    span_max::Dict{String,Float64}
    firstflower_min::Dict{String,Float64}
    firstflower_max::Dict{String,Float64}
    floron::Dict{String,Float64}
    floroff::Dict{String,Float64}
    seednumber_min::Dict{String,Float64}
    seednumber_max::Dict{String,Float64}
    seedon::Dict{String,Float64}
    seedoff::Dict{String,Float64}
    bankduration_min::Dict{String,Float64}
    bankduration_max::Dict{String,Float64}
    b0grow::Dict{String,Float64}
    b0germ::Dict{String,Float64}
    b0mort::Dict{String,Float64}
    temp_opt::Dict{String,Float64}
    temp_tol::Dict{String,Float64}
    fitness::Dict{String, Float64}
end

SppRef() = SppRef(String[],
	          Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict(),
		  Dict())

# Minimal and maximal trait values, which control microevolution, are stored in `traitranges::TraitRanges`
mutable struct TraitRanges
    compartsize::Dict{String,Array{Float64,1}}
    span::Dict{String,Array{Int64,1}}
    firstflower::Dict{String,Array{Int64,1}}
    floron::Dict{String,Array{Int64,1}}
    floroff::Dict{String,Array{Int64,1}}
    seednumber::Dict{String,Array{Int64,1}}
    seedon::Dict{String,Array{Int64,1}}
    seedoff::Dict{String,Array{Int64,1}}
    bankduration::Dict{String,Array{Int64,1}}
end

TraitRanges() = TraitRanges(Dict(),
		            Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict(),
			    Dict())

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
end