command <- function(simID, rseed, cluster, landf, connectsf, sppf, insectsf, disturb, areafile = NULL, duration, tout, tempf, outputat = "EDoutputs", timemsg = "false") {
  
  # generate arguments to be parsed in Julia format 
  parseformat <- function(x, y){
    if(missing(y)){
      paste('\"',ifelse(is.null(x),"",x),'\"', sep = "") #deal with single(not inside folder) or NULL arguments
    }else{
      paste('\"',file.path(y,x),'\"', sep = "")
    }
  }
  
  if (cluster == "hpc") {
    # HPC directories
    outputat <- file.path("/home", outputat)
    EDdir <- file.path("/home/model")
    inputsdir <- file.path("/home/ubuntu/Gaia/inputs")
    Juliadir <- file.path("/home/ubuntu/builds/julia-9d11f62bcb/bin/julia")
  } else {
    # Gaia directories
    EDdir <- file.path("/home/luf74xx/Dokumente/model")
    inputsdir <- file.path(EDdir,"inputs")
    Juliadir <- file.path("/home/luf74xx/builds/julia-d386e40c17/bin/julia")
    outputat <- file.path("EDoutputs")
  }
  command <- paste(
    # julia folder and main
    Juliadir,
    "main.jl", 
    # sim name,
    "--simID",
    parseformat(simID),
    # rseed
    "--rseed",
    parseformat(rseed),
    # outputat
    "--outputat",
    parseformat(outputat),
    "--landconfig",
    # land file
    parseformat(file.path(inputsdir,landf)),
    "--connects",
    # land file
    parseformat(file.path(inputsdir,connectsf)),
    "--spinput", 
    # spp file
    parseformat(file.path(inputsdir,paste(sppf))),
    # insects file
    "--insect",
    parseformat(file.path(inputsdir,insectsf)),
    # type of disturbance
    "--disturb",
    parseformat(disturb),
    ## disturbance file
    "--areafile",
    parseformat(areafile), 
    # sim duration
    "--timesteps",
    parseformat(duration),
    # output frequency
    "--tout",
    parseformat(tout),
    # temperature file
    "--temp_ts",
    parseformat(file.path(inputsdir,tempf)),
    "--timemsg",
    parseformat(timemsg), 
    sep = " ")
  if (cluster == "gaia") {
    scriptname <- paste(simID,".sh", sep = "")
    file.create(file.path(EDdir,scriptname))
    bashscript <- file(file.path(EDdir,scriptname))
    writeLines(c("#!/bin/bash",
                 "#SBATCH -n 12 #number of cores",
                 "#SBATCH --mem-per-cpu=8G",
                 "#SBATCH --output=ED-%j.o",
                 "#SBATCH --error=ED-%j.e",
                 "",
                 command),
               bashscript)
  }else{
    scriptname <- paste(simID,"hpc.txt", sep = "")
    file.create(file.path(file.path("/home/luf74xx/Dokumente/model"),scriptname )) # I create the scripts on Gaia, not in the HPC
    hpcscript <- file(file.path(file.path("/home/luf74xx/Dokumente/model"),scriptname))
    writeLines(command,hpcscript)
    }
}

#command(simID = "testingcommandgen",rseed = 100,cluste="hpc",landf="whatever",sppf="whatever",insectsf="whatever",disturb="whatever",areafile = NULL,duration=10,tout=10,tempf="whatever",outputat = "EDoutputs", timemsg = "false")
