                                        # Subset species from the Göttingen pool
goetspp <- function(simID, rseed, mode, richp = NULL, spplist = NULL){

    library(tidyverse)
    options(scipen=999)

    set.seed(rseed)
    
    EDdir <- file.path("/home/luf74xx/model")
    inputsdir <- file.path(EDdir,"inputs")

                                        # get mean values
    ## table with species traits from leda and Weiss classification
    spptraits <- read.table(file.path(inputsdir, "goettingenEDMtraits.csv"), header = TRUE, sep = ",")
    traits = c()
    
    if (mode == "spplist"){
                                        # read species lists
        spps <- read.table(file.path(inputsdir,spplist), header = TRUE, sep = ",")
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

    ## dispersal kernel
    kernel <- traits$kernel

    ## clonality trait
    clonal <- sample(c("true","false"), richp, replace = TRUE);

    ## get trait mean values and sd of LEDA/Weiss traits
    emass <- traits$emass
    
    elong_mean <- traits$elong_mean
    elong_sd <- traits$elong_sd
    
    floron_mean <- traits$floron_mean
    floron_sd <- traits$floron_sd
    
    floroff_mean <- traits$floroff_mean
    floroff_sd <- traits$floroff_sd

    seedon_mean <- traits$seedon_mean
    seedon_sd <- traits$seedon_sd

    seedoff_mean <- traits$seedoff_mean
    seedoff_sd <- traits$seedoff_sd
    
    span_mean <- traits$span_mean
    span_sd <- traits$span_sd

    ## MTE rates
    b0grow_mean <- traits$b0grow_mean
    b0grow_sd <- traits$b0grow_sd
    
    b0m_mean <- traits$b0m_mean
    b0m_sd <- traits$b0m_sd 

    b0germ_mean <- traits$b0germ_mean
    b0germ_sd <- traits$b0germ_sd
    
    maxseedn_mean <- traits$maxseedn_mean
    maxseedn_sd <- traits$maxseedn_sd

    maxmass <- traits$maxmass


    ## output species id used in EDM
    spEDMid <- data.frame(sp = traits$species, id = paste(rep("p",richp), 1:richp, sep = "-"));
    write.csv(spEDMid, file.path(inputsdir,paste(simID, "ids.csv", sep = "")), row.names = FALSE)

### Write file
    spptable <- data.frame(sp_id = spEDMid$id,
                           kernel = kernel, # R has a function named kernel. Therefore, only the column can take that name
                           clonal = clonal,
                           emass = emass,
                           elong_mean = elong_mean,
                           elong_sd = elong_sd,
                           floron_mean =round(floron_mean,0),
                           floron_sd = round(floron_sd,0),
                           floroff_mean = round(floroff_mean,0),
                           floroff_sd = round(floroff_sd,0),
                           seedon_mean = round(seedon_mean,0),
                           seedon_sd = round(seedon_sd,0),
                           seedoff_mean = round(seedoff_mean,0),
                           seedoff_sd = round(seedoff_sd,0),
                           span_mean = span_mean,
                           span_sd = span_sd,
                           b0grow_mean = b0grow_mean,
                           b0grow_sd = b0grow_sd,
                           b0m_mean = b0m_mean,
                           b0m_sd = b0m_sd,
                           b0germ_mean = b0germ_mean,
                           b0germ_sd = b0germ_sd,
                           maxseedn_mean = maxseedn_mean,
                           maxseedn_sd = maxseedn_sd,
                           maxmass = maxmass)
    
    write.csv(spptable, file.path(inputsdir,paste(simID, ".csv", sep = "")), row.names = FALSE)

    return(spptable)
}

                                        # Creating a random set of species
createsppfile <- function(rseed, richp, sd, simID, b0g, b0em, b0jm, b0am, b0ag, b0jg) {
    set.seed(rseed)
    options(scipen=999)
                                        # sprich = plant species richness
    sp_id <- paste(rep("p",richp), 1:richp, sep = "-");
    
    e_min <- sample(c(0.0001, 0.0003, 0.001), richp, replace = TRUE);
    e_max <- e_min 
    kernels <-sample(c("a","w","wa"), richp, replace = TRUE);
    kernels_sd <- rep(0.1,richp);
    b0g_min <- runif(richp, b0g*0.9,b0g*1.1); #param
    b0g_max <- b0g_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0g*0.9,b0g*1.1)), b0g*0.9,b0g*1.1));
    b0em_min <- runif(richp,b0em*0.9,b0em*1.1);
    b0em_max <- b0em_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0em*0.9,b0em*1.1)), b0em*0.9,b0em*1.1));
    b0jm_min <- runif(richp,b0jm*0.9,b0jm*1.1);
    b0jm_max <- b0jm_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0jm*0.9,b0jm*1.1)), b0jm*0.9,b0jm*1.1));
    b0am_min <- runif(richp, b0am*0.9, b0am*1.1);
    b0am_max <- b0am_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0am*0.9, b0am*1.1)), b0am*0.9, b0am*1.1));
    b0jg_min <- runif(richp, b0jg*0.9, b0jg*1.1); 
    b0jg_max <- b0jg_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0jg*0.9, b0jg*1.1)), b0jg*0.9, b0jg*1.1));
    b0ag_min <- runif(richp, b0ag*0.9, b0ag*1.1);
    b0ag_max <- b0ag_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0ag*0.9, b0ag*1.1)), b0ag*0.9, b0ag*1.1));
    sestra <- sample(c("true","false"), richp, replace = TRUE);
    dyad <- rep(0, richp);
    floron_min <- rep(12, richp)#runif(richp,12,24);
    floron_max <- floron_min + 1 #%>% map_dbl(function(x) ifelse(x < round(mean(c(12,24)),0), 12,24));
    floroff_min <- rep(24, richp)#runif(richp,4,12); #duration
    floroff_max <- floroff_min + 1#%>% map_dbl(function(x) ifelse(x <= mean(c(4,12)), 4,12)); #duration
    seedon_min <- rep(25, richp)#runif(richp,12,24);
    seedon_max <- seedon_min + 1 #%>% map_dbl(function(x) ifelse(x <= mean(c(12,24)), 12,24)); 
    seedoff_min <- rep(37, richp)#runif(richp,4,12); #duration
    seedoff_max <-  seedoff_min + 1#%>% map_dbl(function(x) ifelse(x <= mean(c(4,12)), 4,12)); #duration
    max_mass <- rep(0.0,richp) # calculated once the individual has its see size value
    max_span_min <- rep(1, richp);
    max_span_max <- rep(50, richp);
    first_flower_min <- as.numeric(floor(1.962*max_span_min + 77.24))
    first_flower_max <- as.numeric(floor(1.962*max_span_max + 77.24))
    mass_mu_min <- 0.1*max_mass;
    mass_mu_max <- rep(0.0001, richp);
    abund <-ceiling(runif(richp,100,500));
                                        # Standard deviation valu
    
    spptable <- data.frame(sp_id = sp_id,
                           kernel = kernels, # R has a function named kernel. Therefore, only the column can take that name
                           kernel_sd = kernels_sd,
                           e_mu = e_min,
                           e_sd = e_max,
                           b0g = b0g_min,
                           b0g_sd = b0g_max,
                           b0em = b0em_min,
                           b0em_sd = b0em_max,
                           b0jm = b0jm_min,
                           b0jm_sd = b0jm_max,
                           b0am = b0am_min,
                           b0am_sd = b0am_max,
                           b0jg = b0jg_min,
                           b0jg_sd = b0jg_max,
                           b0ag = b0ag_min,
                           b0ag_sd = b0ag_max,
                           sestra = sestra,
                           dyad = dyad,
                           floron =round(floron_min,0),
                           floron_sd = round(floron_max,0),
                           floroff = round(floroff_min,0),
                           floroff_sd = round(floroff_max,0),
                           seedon = round(seedon_min,0),
                           seedon_sd = round(seedon_max,0),
                           seedoff = round(seedoff_min,0),
                           seedoff_sd = round(seedoff_max,0),
                           max_mass = max_mass,
                           first_flower = first_flower_min,
                           first_flower_sd = first_flower_max,
                           max_span = max_span_min,
                           max_span_sd = max_span_max,
                           mass_mu = mass_mu_min,
                           mass_sd = mass_mu_max,
                           abund = abund)
    
    write.csv(spptable, file.path("/home/luf74xx/Dokumente/model/inputs",paste(simID, ".csv", sep = "")), row.names = FALSE)
    return(spptable)
}
