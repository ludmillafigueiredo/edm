if(Sys.info()["user"] == "ludmilla"){
  EDdir <- file.path("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model")
  EDdocsdir <- file.path(EDdir, "models_docs")
  traitsdir <- file.path("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis","functional_types")
}else if (Sys.info()["user"] == "ubuntu") {
  EDdir <- file.path("/home/ubuntu/model")
  traitsdir <- file.path(EDdir, "inputs")
}else{
  EDdir <- file.path("/home/luf74xx/model")
  EDdocsdir <- file.path("/home/luf74xx/Dokumente/model_docs")
  traitsdir <- file.path(EDdocsdir,"functional_types")
}
inputsdir <- file.path(EDdir, "inputs")
# Set up directory and environement to store analysis
analysEDdir <- file.path(outputsdir, paste(parentsimID, "analysED", sep = "_"))
dir.create(analysEDdir)

if(nreps > 1 | missing(nreps)){
  repfolder <- paste(file.path(outputsdir, parentsimID), 1:nreps, sep = "_"); #parentsimID
}else{
  repfolder <- outputsdir
}
