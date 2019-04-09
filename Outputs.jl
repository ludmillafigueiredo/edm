"""
This module contains functions related to outputting raw results from the model, as well as default derived analysis.
"""
module Outputs

using RCall

export analysED

"""
    analsyED()
Run R script of analysis after the model finishes the simulation
"""
function analysED(settings)
    outdir = settings["outputsat"]
    @rput outdir
    EDdir = pwd()
    @rput EDdir
    nreps = settings["nreps"]
    @rput nreps
    disturbance = settings["disturbtype"]
    @rput disturbance

    # run analysis
    R"analysED.R"
    
end

#end of module
end
