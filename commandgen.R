command <- function(cluster, EDdir, inputsdir, Juliadir, scriptsdir, simID, rseed, outputsdir, spinput, insectsf, initiallandf, disturbtype, landmode, disturbland, tdist, timesteps, tout, tempf, timemsg = "false", nreps = NULL) {
  
  # Generate arguments in Julia parse-format 
  parseformat <- function(x, y){
    if(missing(y)){
      paste('\"',ifelse(is.null(x),"",x),'\"', sep = "") #deal with single(not inside folder) or NULL arguments
    }else{
      paste('\"',file.path(y,x),'\"', sep = "")
    }
  }
  
  # Set up directories
  ## Julia and model folders change according to cluster
  if (cluster == "hpc") {
    # HPC directories
    if (missing(EDdir)){
      EDdir <- file.path("/home/ubuntu/model")
    }
    if (missing(Juliadir)){
      Juliadir <- file.path("/home/ubuntu/build/julia-9d11f62bcb/bin/julia")
    }
  } else {
    # Gaia directories
    if (missing(EDdir)){
      EDdir <- file.path("/home/luf74xx/model")
    }
    
    Juliadir <- file.path("/home/luf74xx/builds/julia-d386e40c17/bin/julia")
  }
  ## input/output directories do not depend on the cluster
  if (missing(inputsdir)){
    inputsdir <- file.path(EDdir, "inputs")
  }
  if (missing(outputsdir)){
    outputsdir <- file.path(EDdir, "EDoutputs")
  }
  ## all scripts are written in Gaia, independently of where they are executed
  if (missing(scriptsdir)){
    scriptsdir <- file.path("/home/luf74xx/model/commandscripts")
  }
  
  # Create command
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
    parseformat(outputsdir),
    "--spinput", 
    # spp file
    parseformat(file.path(inputsdir,paste(spinput))),
    "--insect",
    # insects file
    parseformat(file.path(inputsdir,insectsf)),
    "--initialland",
    # initial raster file
    parseformat(file.path(inputsdir,initiallandf)),
    "--disturbtype",
    # land file
    parseformat(disturbtype),
    # type of disturbance
    "--landmode",
    parseformat(landmode),
    # sim duration
    "--timesteps",
    parseformat(timesteps),
    # output frequency
    "--tout",
    parseformat(tout),
    # temperature file
    "--temp_ts",
    parseformat(file.path(inputsdir,tempf)),
    "--timemsg",
    parseformat(timemsg), 
    sep = " ")
  ## non-obligatory arguments
  if(!is.null(disturbland)){
    command <- paste(command, 
                     "--disturbland",
                     ## disturbance file
                     parseformat(file.path(inputsdir, disturbland)),
                     sep = " ")
  }
  if (!is.null(tdist)){
    command <- paste(command,
                     "--tdist",
                     ## time step of disturbance
                     parseformat(file.path(inputsdir, tdist)),
                     sep = " ")
  }
  if(nreps > 1){
    command <- paste(command,
                     "--nreps",
                     parseformat(nreps),
                     sep = " ")
  }

  # Write command in script for safe-keeping 
  if (cluster == "gaia") {
    scriptname <- paste(simID,".sh", sep = "")
    file.create(file.path(scriptsdir,scriptname))
    bashscript <- file(file.path(scriptsdir,scriptname))
    writeLines(c("#!/bin/bash",
                 "#SBATCH -n 12 #number of cores",
                 "#SBATCH --mem-per-cpu=8G",
                 "#SBATCH --output=ED-%j.o",
                 "#SBATCH --error=ED-%j.e",
                 "",
                 command),
               bashscript)
  }else{
    scriptname <- paste(simID,"hpc.sh", sep = "")
    file.create(file.path(scriptsdir,scriptname)) # I create the scripts on Gaia, not in the HPC
    hpcscript <- file(file.path(file.path(scriptsdir,scriptname)))
    writeLines(c("#!/bin/bash",
                   command,
                   "done"),
                 hpcscript)
  }
  
  
  # bash file to run replicates
  
  #!/bin/bash
  #reps=($(seq 1 29))
  #for i in `seq 1 29`; do
  #simID="v03g74_$i"
  #randomseed="$i"
  #/home/ubuntu/build/julia-9d11f62bcb/bin/julia main.jl --simID $simID --rseed $randomseed --outputat "/home/ubuntu/EDoutputs" --spinput "/home/ubuntu/model/inputs/v03g_66.csv" --insect "/home/ubuntu/model/inputs/insects.csv" --initialland "/home/ubuntu/model/inputs/v03g66initialland.jl" --disturbtype "frag" -landmode "artif" --tdist "/home/ubuntu/model/inputs/tdist52.csv" --timesteps "260" --landbuffer "/home/ubuntu/model/inputs/landbuffer66.jl" --tout "13" --temp <- ts "/home/ubuntu/model/inputs/temp1917_2017.csv" --timemsg "false" --disturbland "/home/ubuntu/model/inputs/v03g_74_disturbland.jl"
  #done
  
}

# Testing lines:

# command(cluster = "gaia",
#         simID = "testingcommandgen", 
#         rseed = 100,
#         spinput = "v03g_59.csv", 
#         insectsf = "insects.csv", 
#         initiallandf = "", 
#         disturbtype = "none", 
#         landmode = "artif", 
#         disturbland = NULL, 
#         tdist = NULL, 
#         timesteps = 5200, 
#         tout = 13, 
#         tempf = "temp1917_2017.csv", 
#         timemsg = "false",  
#         reps = TRUE, 
#         nreps = 30)
# command(cluster = "hpc",
#         simID = "testingcommandgen",
#         rseed = 100,
#         spinput = "v03g_59.csv",
#         insectsf = "insects.csv",
#         initiallandf = "",
#         disturbtype = "none",
#         landmode = "artif",
#         disturbland = NULL,
#         tdist = NULL,
#         timesteps = 5200,
#         tout = 13,
#         timemsg = "false",
#         tempf = "temp1917_2017.csv",
#         reps = TRUE,
#         nreps = 30)
