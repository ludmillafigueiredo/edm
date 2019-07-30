# Subset species from the Göttingen pool

goetspp <- function(inputID, rseed, mode, richp = NULL, spplist = NULL){

    library(tidyverse)
    options(scipen=999)

    set.seed(rseed)
    
    # Define folders
    if(Sys.info()["user"] == "ludmilla"){
      EDdir <- file.path("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/")
      EDdocsdir <- file.path(EDdir, "models_docs/")
      traitsdir <- file.path("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis","functional_types")
    }else{
      EDdir <- file.path("/home/luf74xx/model")
      EDdocsdir <- file.path("/home/luf74xx/Dokumente/model_docs")
      traitsdir <- file.path(EDdocsdir,"functional_types")
    }
    inputsdir <- file.path(EDdir, "inputs")

                                        # get mean values
    ## table with species traits from leda and Weiss classification
    spptraits <- read_csv(file.path(traitsdir, "goetspp_EDMtraits.csv"),
                          col_names = TRUE)
    traits = c()
    
    if (mode == "spplist"){
                                        # read species lists
        spps <- read_table(file.path(spplist), header = TRUE, sep = ",")
                                        # select spp from the list that have known trait values
        traits <- spptraits %>% filter(species %in% spps$species)
        richp <- length(traits$species)
                                        # TODO set up more than one fragment and use map() for more than one
    }else{
                                        # randomly select species for each fragments
                                        # TODO set up more than one fragment and use map() for more than one
        traits <- sample_n(spptraits, richp, replace = FALSE) #richp must have been given 
    }

    ## initial abundances
    abund <-ceiling(runif(richp,20,100))

    ## output species id used in EDM
    spEDMid <- data.frame(sp = traits$sp,
                          id = paste(rep("p", richp), 1:richp, sep = "-"));
    write.csv(spEDMid, file.path(inputsdir, "spp_input", paste(inputID, "ids.csv", sep = "")), row.names = FALSE)

### Write file
    spptable <- data.frame(sp_id = spEDMid$id,
                           abund = abund,
                           select(traits, -sp))
    
    write.csv(spptable, file.path(inputsdir,paste(inputID, "_sppinput.csv", sep = "")), row.names = FALSE)

    return(spptable)
}