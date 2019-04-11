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
function analysED(settings, results_folder)

    parentsimID = settings["simID"]
    @rput parentsimID
    EDdir = pwd()
    @rput EDdir
    nreps = settings["nreps"]
    @rput nreps
    disturbance = settings["disturbtype"]
    @rput disturbance

    outdir = results_folder
    @rput outdir
    
    # run analysis
    R"analysED.R"
    
end

#end of module
end
