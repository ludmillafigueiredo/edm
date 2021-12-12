#!/Usr/bin/env julia

const EDMdir = abspath(joinpath(pwd(), "src")) # absolute path because it is sent to R - avoid mistakes
push!(LOAD_PATH, EDMdir)
#cd(EDDir)

include(joinpath(EDMdir,"setup_julia.jl")) # for "using" statements
include(joinpath(EDMdir,"constants.jl"))
include(joinpath(EDMdir,"entities.jl"))
include(joinpath(EDMdir,"auxiliary.jl"))
include(joinpath(EDMdir,"inputs.jl"))
include(joinpath(EDMdir,"initialisation.jl"))
include(joinpath(EDMdir,"submodels_lifecycle.jl"))
include(joinpath(EDMdir,"submodels_disturbances.jl"))
include(joinpath(EDMdir,"outputs.jl"))
include(joinpath(EDMdir,"sim_log.jl"))
include(joinpath(EDMdir,"unit_tests.jl"))
include(joinpath(EDMdir,"scheduling.jl"))

start_time = now()

# run simulation
run_scheduling(settings, management_counter, land_pars, poll_pars, K, T, mean_annual, plants, landscape)

end_time = now()

println("Simulation time: ", end_time-start_time)

start_time = now()

"""
# analyse results
analysED(settings, land_pars, poll_pars)

end_time = now()

println("Analyze time: ", end_time-start_time)
"""
