#### Functions ####

#' Organize individual-based output into csv files of state variables (only available for juveniles and adults).
#' @param simID The simulation ID
#' @param outputsdir The path to the outputs folder, if not default
#' @param EDMdir The path to EDM, if not default
format_statevars <- function(outputsdir, simID, analysEDMdir){

  spp_seedmass <- read_csv(file.path(outputsdir, "sppref_traitvalues.csv"), 
                           col_names=TRUE) %>%
		  dplyr::select(sp, seedmass) 

  raw <- read_tsv(file.path(outputsdir, "statevars_ind.txt"), col_names = TRUE,
                       col_types = cols(week = col_integer(),
                                        id = col_character(),
                                        stage = col_character(),
                                        location = col_character(),
                                        sp = col_character(),
                                        kernel = col_character(),
                                        clonality = col_logical(),
                                        pollen_vector = col_character(),
                                        self_failoutcross = col_character(),
                                        self_proba = col_double(),
                                        compartsize = col_double(),
                                        span = col_integer(),
                                        firstflower = col_integer(),
                                        floron = col_integer(),
                                        floroff = col_integer(),
                                        seednumber = col_integer(),
                                        seedon = col_integer(),
                                        seedoff = col_integer(),
                                        bankduration = col_integer(),
                                        age = col_integer(),
                                        leaves = col_double(),
                                        stem = col_double(),
                                        root = col_double(),
                                        repr = col_double(),
                                        mated = col_logical()));
    # clean it
    ## take parentheses out of location column ("()")
    loc <- gsub("[\\(|\\)]", "", raw$location)
    loc <- as.data.frame(matrix(unlist(str_split(loc, ",")), ncol = 2, byrow = T))
    names(loc) = c("xloc", "yloc")
    ## complete and clean table
    outdata <- bind_cols(dplyr::select(raw, -location), loc)
    rm(loc)
    ## add colum containing values of seedmass
    ## get rid of 'mated': reproduction is better assessed in "pollination output.
    ## add column identifying replicate
    statevars <- inner_join(outdata, spp_seedmass)%>%
    	       select(-mated)%>%
	       mutate(simID = simID)
      
    ## write single files in its folders
    write_csv(statevars,
	      file.path(analysEDMdir, "statevars.csv"))
  
  return(statevars)
  
}

#' Format seed output and id simulation it belongs to.
#' @param simID simulation ID
#' @outputsdir directory where simulation outputs were writen.
#' @analysEDMdir directory where analysis are writen.
format_seedprod <- function(outputsdir, simID, analysEDMdir){
  
  seed_production <- read_tsv(file.path(outputsdir, "offspringproduction.csv"),
    		    	     col_names = TRUE,
			     col_types = cols(week = col_integer(),
                                        stage = col_character(),
                                        sp = col_character(),
                                        abundance = col_double(),
                                        mode = col_character()))%>%
      mutate(simID = simID)
    
  ## write single files in its folders
  write_csv(seed_production,
	    file.path(analysEDMdir, "seedproduction.csv"))

  return(seed_production)

}

#' Extract and plot species-specific mean and sd
#' of biomass compartments of juveniles and adults.
#' @param outdatalist
#' @param plotit Boolean specifying whether graph shouuld be plotted or not
biomass_allocation <- function(statevars, analysEDMdir, stages = factor(c("a", "j", "s"))){
  
  ## get biomasses allocated to each compartment
  allocation_summary <- statevars %>%
    dplyr::select(simID, week, id, stage, seedmass, leaves, stem, root, repr)%>%
    filter(stage %in% c("j", "a"))%>%
    gather(key = compart, value = mass, leaves:repr) %>%
    group_by(simID, week, stage, seedmass, compart) %>%
    summarize(mass_mean = mean(mass),
              mass_sd = sd(mass)) %>%
    ungroup()

  allocation_plot <- allocation_summary %>%
    ggplot(aes(x = week, y = mass_mean, colour = compart))+
    geom_line()+
    geom_errorbar(aes(ymin = mass_mean - mass_sd,
    		      ymax = mass_mean + mass_sd),
		      width = 0.1)+
    facet_grid(stage ~ seedmass, 
               scales = "free_y")
  
  write_csv(allocation_summary,
	      file.path(analysEDMdir, "allocation_summary.csv"))    
  
  return(allocation_plot)
}

#' Plot growth curve
growthcurves <- function(statevars, analysEDMdir){
    
  growthcurve_tab <- statevars%>%
    dplyr::select(week, id, stage, seedmass, simID, leaves, stem, root, repr)%>%
    filter(stage %in% c("j", "a"))%>%
    mutate(total = leaves + stem + root + repr) %>%
    group_by(week, stage, seedmass, simID) %>%
    summarize(totalmass_mean = mean(total),
              totalmass_sd = sd(total)) %>%
    ungroup() 

  growthcurve_plot <- growthcurve_tab%>%
    mutate(seedmass = as.factor(seedmass)) %>%
    ggplot(aes(x = week, y = totalmass_mean, group = seedmass, color = seedmass))+
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin = totalmass_mean - totalmass_sd,
    		      ymax = totalmass_mean + totalmass_sd),
		      width = 0.1)+
    facet_wrap(~stage, ncol = 1,
               scales = "free_y")+
    labs(title = "Individual growth curves")+
    theme(legend.position = "bottom")

    write_csv(growthcurve_tab,
	      file.path(analysEDMdir, "growthcurve.csv"))    
    
    return(growthcurve_plot)
}

#' Extract and plot population abundances
pop_abundances <- function(statevars, analysEDMdir){
  
  ## extract relevant variables
  pop_vars <- statevars%>%
    dplyr::select(week, id, sp, stage, simID)%>%
    ungroup()
  
  ## summarize abundaces/species/week
  spabund_plot <- pop_vars%>%
    group_by(week, sp, simID)%>%
    summarize(abundance = n())%>%
    ungroup()%>%
    ggplot(aes(x = week, y = abundance,
    	       group = sp, color = sp))+
    geom_line()+
    geom_point()+
    labs(x = "Year", y = "Abundance",
         title = "Species abundance variation")+
    theme(legend.position = "none")
    
  write_csv(pop_vars,
	    file.path(analysEDMdir, "pop_variables.csv"))
  
  return(list(tab = pop_vars, plot = spabund_plot))
}

#' Extract and plot populations structure
pop_structure <- function(pop_vars, seed_production, simID, analysEDMdir){
  
  ## Population strucuture: absolute abundances
    absltstruct_tab <- pop_vars%>%
      group_by(week, sp, stage)%>%
      summarize(abundance = n())%>%
      ungroup()%>%
      bind_rows(., dplyr::select(seed_production, -mode))%>% #merge offspring file
      ungroup()
    
    ## Population strucuture: relative abundances
    rltvstruct_tab <- pop_vars%>%
      group_by(week, sp, stage)%>%
      summarize(abundance = n())%>%
      complete(week, stage, fill = list(abundance = 0))%>%
      group_by(week, sp)%>%
      mutate(total_abund = sum(abundance),
             rltv_abund = abundance/total_abund)%>%
      complete(week, stage, fill = list(rel_abund = 0))%>%
      ungroup()
  
  ## Graphs
  absltstruct_plot <- ggplot(absltstruct_tab,
                             aes(x = week, y = abundance, color = factor(stage)))+
    geom_line()+
    geom_point(size = 1)+
    labs(x = "Week", y = "Abundance",
         title = "Population structure")+
    facet_wrap(~sp, ncol = 8, scale = "free_y")+
    theme(legend.position = "bottom")
  
  rltvstruct_plot <- ggplot(rltvstruct_tab,
                            aes(x = week, y = rltv_abund, color = as.factor(stage)))+
    geom_line()+
    geom_point(size = 1)+
    labs(x = "Week", y = "Relative abundace",
         title = "Population structure")+
    facet_wrap(~sp, ncol = 8)+
    theme(legend.position = "bottom")
    
    
    write_csv(absltstruct_tab,
	    file.path(analysEDMdir, "absolute_structure.csv"))
    write_csv(rltvstruct_tab,
	    file.path(analysEDMdir, "relative_structure.csv"))    

    return(list(absplot = absltstruct_plot, rltvplot = rltvstruct_plot))
  
}

#' Calculate mean species richness per group of size
spp_richness <- function(statevars, pop_vars, disturbance, tdist, analysEDMdir){
  
  ## extract richness from output
  spprichness_tab <- pop_vars %>%
    group_by(week, simID) %>%
    summarize(richness = length(unique(sp))) %>%
    ungroup()
  
  ## create plots
  ## identify when disturbance happens, if it does
  spprichness_plot <- ggplot(spprichness_tab, aes(x = week, y = richness))+
    geom_line()+
    geom_point()+
    labs(x = "Time", y = "Species richness")
    
  if (disturbance != "none"){
      spprichness_plot <- spprichness_plot +
                      geom_vline(xintercept = tdist, linetype = 2, color = "red", alpha = 0.5)      
  }
  
  ## richness per group of size
  groupspprichness_tab  <- statevars%>%
    dplyr::select(week, sp, seedmass, simID)%>%
    ungroup()%>%
    group_by(week, seedmass, simID)%>%
    summarize(richness = length(unique(sp)))%>%
    ungroup()
  
  grouprichness_plot <- ggplot(data = groupspprichness_tab, aes(x = week, y = richness))+
    facet_wrap(~seedmass, scales = "free_y", ncol = 1,
               labeller = labeller(seedmass = size_labels))+
    geom_line()+
    geom_point(size = 1.25)+
    labs(x = "Time", y = "Species richness")
  
  write_csv(spprichness_tab,
	    file.path(analysEDMdir, "spprichness.csv"))

  return(grouprichness_plot)
}

#' Calculate species relative abundance and plot species rank-abundance at a given time-step
rank_abundances <- function(pop_vars, timesteps, analysEDMdir){
  
  spprltvabund_tab <- pop_vars%>%
    filter(week %in% timesteps)%>%
    group_by(week, sp)%>%
    summarize(abundance = n())%>%
    group_by(week)%>%
    mutate(total_abund = sum(abundance),
           rltv_abund = abundance/total_abund)
  
  plot_rankabund <- function(timestep){
    ggplot(filter(spprltvabund_tab, week %in% timestep),
           aes(x = reorder(sp, -rltv_abund), y = rltv_abund))+
      geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 80, size = 10, vjust = 0.5))
  }
  
  rankabund_plots <- timesteps%>%
    map(. %>% plot_rankabund)
  
  rankabund_plot <- plot_grid(plotlist = rankabund_plots, 
                              labels = paste("week", timesteps, sep = " "))
  
  write_csv(spprltvabund_tab,
	    file.path(analysEDMdir, "sppabundances.csv"))
  
  return(rankabund_plot)

}

#' Calculate biomass production
biomass_production <- function(statevars, analysEDMdir){
  
  production_tab <- statevars%>%
    dplyr::select(week, leaves, stem, root, repr, simID)%>%
    mutate(production = (leaves+stem+root+repr)/(10^3))%>%
    ungroup()%>%
    group_by(week, simID)%>%
    summarize(prod_mean = mean(production),
              prod_sd = sd(production))%>%
    ungroup()
  
  production_plot <- ggplot(production_tab, aes(x= week, y = prod_mean))+
    geom_errorbar(aes(ymin = prod_mean - prod_sd,
                      ymax = prod_mean + prod_sd),
                  stat = "identity",
		      width = 0.1)+
    geom_point()+
    geom_line()+
    labs(x = "Week",
         y = "Biomass production (kg)")

  write_csv(production_tab,
	    file.path(analysEDMdir, "biomass_production.csv"))
  
  return(production_plot)

}


#' Life-history events
lifehistory <- function(outputsdir, simID, analysEDMdir){
  
  lifeevents_tab <- read_tsv(file.path(outputsdir, "eventslog.txt"), col_names = TRUE) %>%
  		    mutate(simID = simID)
  lifeevents_plot <- lifeevents_tab%>%
    dplyr::select(-age)%>%
    group_by(week, event, stage, simID)%>%
    summarize(total = n())%>%
    ggplot(aes(x = week, y = total, color = event))+
    geom_line()+
    facet_wrap(~stage, ncol = 1, scales = "free_y")
  
  #metabolicrates_tab <- read_tsv(file.path(outputsdir, "metaboliclog.txt"), col_names = TRUE) %>%
  #  dplyr::select(-age)%>%
  #  gather(c(rate, probability), key = "metric", value = value)%>%
  #  group_by(stage, event, metric)%>%
  #  summarize(mean = mean(value),
  #            sd = sd(value)) %>%
  #  mutate(simID = simID)
  
  write_csv(lifeevents_tab,
	    file.path(analysEDMdir, "lifeevents.csv"))
  
  return(lifeevents_plot)

}

#' Age distribution of stages
agetraits <- function(statevars, analysEDMdir){
  
  meanage_tab <- statevars%>%
    group_by(week, sp, stage, simID)%>%
    summarize(mean_age = mean(age),
              sd_age = sd(age),
              mean_span = mean(span),
              sd_span = sd(span),
              mean_firstflower = mean(firstflower),
              sd_firstflower = sd(firstflower))%>%
    ungroup()
  
  meanage_juvplot <- meanage_tab%>%
    filter(stage == "j")%>%
    ggplot(aes(x=week, y=mean_age, colour=stage))+
    geom_hline(aes(yintercept = mean_firstflower), colour = "gold", alpha = 0.1)+
    geom_vline(aes(xintercept = mean_firstflower), colour = "gold", alpha = 0.1)+
    geom_errorbar(aes(ymin=mean_age-sd_age, ymax=mean_age+sd_age),
		      width = 0.1)+
    geom_line()
    
  meanage_adtplot <- meanage_tab%>%
    filter(stage == "a")%>%
    ggplot(aes(x=week, y=mean_age, colour=stage))+
    geom_hline(aes(yintercept = mean_span), colour = "grey", alpha = 0.1)+
    geom_vline(aes(xintercept = mean_span), colour = "grey", alpha = 0.1)+
    geom_errorbar(aes(ymin=mean_age-sd_age, ymax=mean_age+sd_age),
		      width = 0.1)+
    geom_line()
  
  meanage_plot <- plot_grid(meanage_juvplot, meanage_adtplot,
                            ncol = 2, nrow = 1)

  return(list(tab = meanage_tab, jplot = meanage_juvplot, aplot = meanage_adtplot))

}

#' Seed production and seed bank
#seeddynamics <- function(statevars, seed_production, analysEDMdir){
#  
#  seeddyn_tab <- statevars%>%
#    filter(stage == "s")%>%
#    dplyr::select(week, sp, stage, simID)%>%
#    group_by(week, sp, stage, simID)%>%
#    summarize(abundance = n())%>%
#    ungroup()%>%
#    mutate(stage = "seedbank")%>%
#    bind_rows(dplyr::select(seed_production, -mode))
#  
#  seeddyn_plot <- ggplot(seeddyn_tab, aes(x = week, y = abundance))+
#    geom_line(aes(group = sp, colour = stage))+
#    geom_point()
#    
#  return(list(a = seeddyn_tab, b = seeddyn_plot))
#}

