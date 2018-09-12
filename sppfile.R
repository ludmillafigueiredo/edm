# Creating a random set of species
createsppfile <- function(rseed, richp, sd, simID, tdist, b0g, b0em, b0am, b0ag, b0jg) {
  set.seed(rseed)
  # sprich = plant species richness
  sp_id <- paste(rep("p",richp), 1:richp, sep = "-");
  
  if (tdist == "normal") {
    # Mean trait values: LEDA, parameterization or rdmly assigned
    e_mu <- runif(richp, min = 0.001, max = 0.07);
    e_sd <- sd * e_mu
    kernels <-sample(c("a","w","wa"), richp, replace = TRUE);
    kernels_sd <- rep(sd,richp)
    b0g <- rep(1502900000000, richp); #param
    b0g_sd <- sd * b0g
    b0em <- rep(1502900000000, richp);
    b0em_sd <- sd * b0em
    b0am <- rep(47527000000, richp);
    b0am_sd <- sd * b0am
    b0jg <- rep(67026000000, richp);
    b0jg_sd <- sd * b0jg
    b0ag <- rep(335130000, richp);
    b0ag_sd <- sd * b0ag
    sestra <- rep(0, richp);
    dyad <- rep(0, richp);
    floron <- rep(12, richp);
    floron_sd <- sd * floron
    floroff <- ceiling(runif(richp, min = 1, max = 12));
    floroff_sd <- sd * floroff
    seedon <- rep(32, richp);
    seedon_sd <- sd * seedon
    seedoff <- ceiling(runif(richp, min = 1, max = 5));
    seedoff_sd <- seedoff
    max_mass <- ((e_mu*10^(3))/2.14)^2
    max_span <- ceiling(runif(richp, min = 1, max = 50));
    max_span_sd <- max_span
    mass_mu <- 0.8*max_mass;
    mass_sd <- rep(0.0001, richp);
    abund <- ceiling(runif(richp,1,20));
  } else {
    e_min <- runif(richp, 0.007, 0.07);
    e_max <- e_min %>% map_dbl(function(x) ifelse(x <= mean(c(0.007, 0.07)), 0.007, 0.07));
    kernels <-sample(c("a","w","wa"), richp, replace = TRUE);
    kernels_sd <- rep(0.1,richp);
    b0g_min <- runif(richp, b0g*0.8,b0g*1.2); #param
    b0g_max <- b0g_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0g*0.8,b0g*1.2)), b0g*0.8,b0g*1.2));
    b0em_min <- runif(richp,b0em*0.8,b0em*1.2);
    b0em_max <- b0em_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0em*0.8,b0em*1.2)), b0em*0.8,b0em*1.2));
    b0am_min <- runif(richp, b0am*0.8, b0am*1.2);
    b0am_max <- b0am_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0am*0.8, b0am*1.2)), b0am*0.8, b0am*1.2));
    b0jg_min <- runif(richp, b0jg*0.8, b0jg*1.2); 
    b0jg_max <- b0jg_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0jg*0.8, b0jg*1.2)), b0jg*0.8, b0jg*1.2));
    b0ag_min <- runif(richp, b0ag*0.8, b0ag*1.2);
    b0ag_max <- b0ag_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0ag*0.8, b0ag*1.2)), b0ag*0.8, b0ag*1.2));
    sestra <- rep(0, richp);
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
    abund <-ceiling(runif(richp,1,20));
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
    
  }
  write.csv(spptable, file.path("/home/luf74xx/Dokumente/model/inputs",paste(simID, ".csv", sep = "")), row.names = FALSE)
  return(spptable)
}

# Subset species from the Göttingen pool
goetspp <- function(rseed,richp,mode,simID){
  set.seed(rseed)
  
  EDdir <- file.path("/home/luf74xx/Dokumente/model")
  inputsdir <- file.path(EDdir,"inputs")
  
  # Get *seed size*, *first age of flowering* and *maximum life span* from LEDA
  # fill the ones that are NA for span or first age
  require(tidyverse)
  spp <- read.table(file.path(inputsdir,"traitsleda.csv"), header = TRUE, sep = ",")%>%
    select(sp,seedm,kernels,span,firstflower)%>%
    mutate(seedm = 0.001*seedm)%>%
    mutate(kernels = case_when(
      kernels == "wind" ~ "w",
      kernels == "ant" ~ "a",
      kernels == "windant" ~ "wa"
    ))%>%
    mutate(span = case_when(
      is.na(span) ~ if_else(firstflower < 75.0, as.integer(1), as.integer(ceiling((firstflower - 77.24)/1.962))),
      TRUE ~ span))%>%
    mutate(firstflower = case_when(
        is.na(firstflower) ~ as.numeric(floor(1.962*span + 77.24)),
        TRUE ~ firstflower))%>%
    filter(!is.na(span))#filter those out, for now
    
  if (mode == "frags"){
    # read fragments composition and create species lists
  }else{
    # randomly select species for each fragments
    # TODO set up more than one fragment and use map() for more than one
    traits <- sample_n(spp, richp, replace = FALSE)
  }
    spEDid <- data.frame(sp = traits$sp, id = paste(rep("p",richp), 1:richp, sep = "-"));
    
    # get other traits
    e_min <- 0.8*traits$seedm
    e_max <- 1.2*traits$seedm
    kernels <- traits$kernels
    kernels_sd <- rep(0.1,richp); # doesnt do anything yet
    b0g_min <- runif(richp, b0g*0.8,b0g*1.2); #param
    b0g_max <- b0g_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0g*0.8,b0g*1.2)), b0g*0.8, b0g*1.2));
    b0em_min <- runif(richp,b0em*0.8,b0em*1.2);
    b0em_max <- b0em_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0em*0.8,b0em*1.2)), b0em*0.8, b0em*1.2));
    b0am_min <- runif(richp, b0am*0.8, b0am*1.2);
    b0am_max <- b0am_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0am*0.8, b0am*1.2)), b0am*0.8, b0am*1.2));
    b0jg_min <- runif(richp, b0jg*0.8, b0jg*1.2); 
    b0jg_max <- b0jg_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0jg*0.8, b0jg*1.2)), b0jg*0.8, b0jg*1.2));
    b0ag_min <- runif(richp, b0ag*0.8, b0ag*1.2);
    b0ag_max <- b0ag_min %>% map_dbl(function(x) ifelse(x <= mean(c(b0ag*0.8, b0ag*1.2)), b0ag*0.8, b0ag*1.2));
    sestra <- rep(0, richp);
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
    first_flower_min <- ceiling(0.8*traits$firstflower)
    first_flower_max <- ceiling(1.2*traits$firstflower)
    max_span_min <- ceiling(0.8*traits$span)
    max_span_max <- ceiling(1.2*traits$span)
    mass_mu_min <- 0.1*max_mass;
    mass_mu_max <- rep(0.0001, richp);
    abund <-ceiling(runif(richp,1,20))
    
    spptable <- data.frame(sp_id = spEDid$id,
                           kernel = kernels, # R has a function named kernel. Therefore, only the column can take that name
                           kernel_sd = kernels_sd,
                           e_mu = e_min,
                           e_sd = e_max,
                           b0g = b0g_min,
                           b0g_sd = b0g_max,
                           b0em = b0em_min,
                           b0em_sd = b0em_max,
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
    write.csv(spEDid, file.path("/home/luf74xx/Dokumente/model/inputs",paste(simID, "ids.csv", sep = "")), row.names = FALSE)
    return(spptable)
}

