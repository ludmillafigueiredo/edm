#!/Usr/bin/env julia

const EDMdir = joinpath(pwd(), "src") # TODO test without pwd()
push!(LOAD_PATH, EDMdir)
#cd(EDDir)

include("constants.jl")
include("entities.jl")
include("auxiliary.jl")
include("inputs.jl")
include("initialisation.jl")
include("submodels_lifecycle.jl")
include("submodels_disturbance.jl")
include("outputs.jl")
include("unit_tests.jl")
include("scheduling.jl")

# run simulation
results_folder = run_scheduling(settings, management_counter, landpars, poll_pars, K, T, mean_annual, plants, landscape)

# analyse results
analysED(settings, results_folder)
