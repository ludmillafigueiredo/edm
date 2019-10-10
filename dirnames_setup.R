if(Sys.info()["user"] == "ludmilla"){
  EDdir <- file.path("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/src")
  EDdocsdir <- file.path(EDdir, "models_docs")
  traitsdir <- file.path("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis","functional_types")
}else if (Sys.info()["user"] == "ubuntu") {
  EDdir <- file.path("/home/ubuntu/model/src")
  traitsdir <- file.path(EDdir, "inputs")
}else{
  EDdir <- file.path("/home/luf74xx/model/src")
  EDdocsdir <- file.path("/home/luf74xx/Dokumente/model_docs")
  traitsdir <- file.path(EDdocsdir,"functional_types")
}
inputsdir <- file.path(EDdir, "inputs")

if(nreps > 1 | missing(nreps)){
  repfolder <- paste(file.path(outputsdir, parentsimID), 1:nreps, sep = "_"); #parentsimID
}else{
  repfolder <- outputsdir
}
