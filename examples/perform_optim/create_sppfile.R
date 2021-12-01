# Subset species from the Göttingen pool

create_sppfile <- function(inputID, rseed, mode, richp = NULL, traitsmode, traitsdir, spplist = NULL, inputsdir){

    require(tidyverse)
    options(scipen=999)

    set.seed(rseed)
    
    # Define folders
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
    inputsdir <- file.path("examples/perform_optim")

                                        # get mean values
    ## table with species traits from leda and Weiss classification
    if (traitsmode == "normal"){
      spptraits <- read_csv(file.path(EDdocsdir, "goetspp_EDMtraits_normal.csv"),
                            col_names = TRUE)
    }else{
      spptraits <- read_csv(file.path("docs/TRACE/spp_traits/EDM_spptraits.csv"),
                            col_names = TRUE)
      }
      
    end
    spptraits$species <- str_replace(spptraits$species, " ", "_")
    
    traits = c()
    
    if (mode == "spplist"){
                                        # read species lists
        spps <- read_csv(file.path(EDdocsdir, spplist), col_names = TRUE)
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
    sppinput <- traits%>%
                add_column(abund = ceiling(runif(richp,20,100)),
                           sp_id = paste(rep("p", richp), 1:richp, sep = "-"))

    ## output species id used in EDM            
    write.csv(dplyr::select(sppinput, species, sp_id),
              file.path(inputsdir, paste(inputID, "ids.csv", sep = "")), row.names = FALSE)

### Write file
    write.csv(dplyr::select(sppinput, -species),
              file.path(inputsdir,paste(inputID, "_sppinput.csv", sep = "")), row.names = FALSE)

    return(sppinput)
}
