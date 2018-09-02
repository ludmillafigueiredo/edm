command <- function(simID,rseed,cluster,landf,sppf,insectsf,disturbf,duration,tout,tempf) {
  if (cluster == "hpc") {
    EDdir <- file.path("/home/ubuntu/model")
    Juliadir <- file.path("/home/ubuntu/builds/julia-9d11f62bcb/bin/julia")
  } else {
    EDdir <- file.path("/home/luf74xx/Dokumente/model")
    Juliadir <- file.path("/home/luf74xx/builds/julia-d386e40c17/bin/julia")
  }
  command <- paste(
    # julia folder and main
    Juliadir,
    "main.jl", 
    # arg 1
    "--simID",
    # sim name,
    paste('\"',simID,'\"', sep = ""),
    "--rseed",
    paste('\"',rseed,'\"', sep = ""),
    # arg 2
    "--landconfig",
    # land file
    paste('\"',file.path(EDdir,"inputs",landf),'\"', sep = ""),
    # arg 3
    "--spinput", 
    # spp file
    paste('\"',file.path(EDdir,"inputs",sppf),'\"', sep = ""),
    # arg 4
    "--insect",
    # insects file
    paste('\"',file.path(EDdir,"inputs",insectsf),'\"', sep = ""), 
    # arg 5
    "--disturb",
    # disturbance file
    paste('\"',disturbf,'\"', sep = ""),
    #arg 6
    "--timesteps", 
    # sim duration
    paste('\"',duration,'\"', sep = ""),
    # arg 8
    "--tout",
    # output frequency
    paste('\"',tout,'\"', sep = ""),
    "--temp_ts",
    # arg 9
    # temperature file
    paste('\"',file.path(EDdir,"inputs",tempf),'\"', sep = ""), 
    sep = " ")
  if (cluster == "gaia") {
    file.create(file.path(EDdir,paste(simID,".sh", sep = "")))
    bashscript <- file(file.path(EDdir,paste(simID,".sh", sep = "")))
    writeLines(c("#!/bin/bash",
                 "#SBATCH -n 12 #number of cores",
                 "#SBATCH --mem-per-cpu=8G",
                 "#SBATCH --output=ED-%j.o",
                 "#SBATCH --error=ED-%j.e",
                 "",
                 command),
               bashscript)
  }else{
    return (command)
    }
}