## This is an R script intended at fast analysis. Details on the RNotebook file of same name
##

#### Load packages ####
library(tidyverse);
library(viridis);
library(grid);
library(gridExtra);
library(FactoMineR);
library(factoextra); #might require mvtnorm 1.0-6, which does not require R 3.5
		     #devtools::install_version("mvtnorm", version = "1.0-6",
		     #                          repos = "http://cran.us.r-project.org")
library(corrplot);
library(gganimate); #sudo apt-get install cargo install.packages("gifski")
library(cowplot);

#### Directories and environment to store analysis ####
#source("src/dirnames_setup.R")
analysEDdir <- file.path(outputsdir, paste(parentsimID, "analysED", sep = "_"))
dir.create(analysEDdir)

#### Theme, colours, labels ####
source("src/theme_edm.R")
theme_set(theme_edm())
size_labels = c('0.0001' = "Small",
                '0.0003' = "Medium",
                '0.001' = "Big")
options(scipen = 999)


#### Functions ####

#' Organize individual-based output
#' 
#' @param parentsimID The simulation ID
#' @param nreps Number of replicates
#' @param outputsdir The path to the outputs folder, if not default
#' @param EDdir The path to EDM, if not default
getoutput <- function(parentsimID, repfolder, nreps, outputsdir, EDdir = file.path("~/model")){

  spp_seedmass <- read_csv(file.path(outputsdir, "sppref_traitvalues.csv"), 
                           col_names=TRUE) %>% dplyr::select(sp, seedmass) 

  ## initiliaze list to contain outputs 
  ##outdatalist <- list()
  for(repli in repfolder){
    
    ## get raw outputs
    outraw <- read_tsv(file.path(repli, "statevars_ind.txt"), col_names = TRUE,
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
    loc <- gsub("[\\(|\\)]", "", outraw$location)
    loc <- as.data.frame(matrix(unlist(str_split(loc, ",")), ncol = 2, byrow = T))
    names(loc) = c("xloc", "yloc")
    ## complete and clean table
    outdata <- bind_cols(dplyr::select(outraw, -location), loc)
    rm(loc)
    ## add colum containing values of seedmass
    outdata <- inner_join(outdata, spp_seedmass)
      
    ## write single files in its folders
    write.csv(outdata,
	      file.path(repli, paste(parentsimID, "juvads_statevars.csv", sep = "_")), row.names = FALSE)
    
  }
  
  return(outdata)
  
}

#' Assemble output files of replicates and identify them.
#' Output on juvenile and adult individuals are kept separate because the offspring
#' contains weekly abundances of each species, whereas the juv/adult output is individual-based.
#' @param parentsimID simulation ID
#' @param nreps number of replicates, which identify the folders containing results
orgareplicates <- function(parentsimID, repfolder, nreps){
  
  juvads_allreps_tab <- data.frame()
  seeds_allreps_tab <- data.frame()
  
  for(i in 1:nreps) {
    sim <- paste(parentsimID, i, sep = "_")
    folder <- repfolder[i]
    ## read outputs of juv/adults and of offspring
    outdata_sim <- read_csv(file.path(folder, paste(parentsimID, "juvads_statevars.csv", sep = "_")),
    	       	           col_names = TRUE)%>%
        mutate(repli = as.factor(rep(sim, nrow(.))))
    offspring_sim <- read_tsv(file.path(folder, "offspringproduction.csv"),
    		    	     col_names = TRUE,
			     col_types = cols(week = col_integer(),
                                        stage = col_character(),
                                        sp = col_character(),
                                        abundance = col_double(),
                                        mode = col_character()))%>%
      mutate(repli = as.factor(rep(sim, nrow(.))))
    
    ## fill in cleanoutput object (necessary step because of bind_rows)
    if(length(juvads_allreps_tab) == 0){
      juvads_allreps_tab <- outdata_sim
      seeds_allreps_tab <- offspring_sim
    }else{
      ## identify the lines corresponding to it with its replicate ID
      juvads_allreps_tab <- bind_rows(juvads_allreps_tab, outdata_sim, .id = "repli")
      seeds_allreps_tab <- bind_rows(seeds_allreps_tab, offspring_sim, .id = "repli")
    }
  }
    
  return(list(a = juvads_allreps_tab, b = seeds_allreps_tab))
}

#' Extract and plot species-specific mean and sd
#' of biomass compartments of juveniles and adults.
#' @param outdatalist
#' @param plotit Boolean specifying whether graph shouuld be plotted or not
biomass_allocation <- function(juvads_allreps_tab, stages = factor(c("a", "j", "s"))){
  
  ## get biomasses allocated to each compartment
  # summary within replicates
  mass_summaryrepli <- juvads_allreps_tab%>%
    dplyr::select(week, id, stage, sp, leaves, stem, root, repr, repli)%>%
    filter(stage %in% c("j", "a"))%>%
    group_by(week, stage, sp, repli)%>%
    summarize(mean_reprmass = mean(repr),
              mean_leavesmass = mean(leaves),
              mean_stemmass = mean(stem),
              mean_rootmass = mean(root),
              sd_reprmass = sd(repr),
              sd_leavesmass = sd(leaves),
              sd_stemmass = sd(stem),
              rootmass = sd(root))%>%
    ungroup()
  # summary through replicates
  mass_tab <- juvads_allreps_tab%>%
    dplyr::select(week, id, stage, sp, leaves, stem, root, repr)%>%
    filter(stage %in% c("j", "a"))%>%
    group_by(week, stage, sp)%>%
    summarize(mean_repr = mean(repr),
              mean_leaves = mean(leaves),
              mean_stem = mean(stem),
              mean_root = mean(root),
              sd_repr = sd(repr),
              sd_leaves = sd(leaves),
              sd_stem = sd(stem),
              sd_root = sd(root))%>%
    ungroup()
  
  # Plot biomass in compartments
  # TODO: Factorize it, but beware of the compartment specific column names.
  # ends_with() does not seem to work inside aes
  leaves_plot <- mass_tab%>%
    dplyr::select(week, stage, sp, ends_with("leaves"))%>%
    ggplot(aes(x = week, y = mean_leaves, group = factor(sp), colour = factor(sp)))+
    geom_errorbar(aes(min = mean_leaves - sd_leaves,
                      max = mean_leaves + sd_leaves),
                  stat = "identity")+
    geom_line() +
    geom_point()+
    facet_wrap(~stage, ncol = 2)+
    labs(title = "Biomass allocated to leaves")+
    theme(legend.position = "none")
  
  stem_plot <- mass_tab%>%
    dplyr::select(week, stage, sp, ends_with("stem"))%>%
    ggplot(aes(x = week, y = mean_stem, group = factor(sp), colour = factor(sp)))+
    geom_errorbar(aes(min = mean_stem - sd_stem,
                      max = mean_stem + sd_stem),
                  stat = "identity")+
    geom_line() +
    geom_point()+
    facet_wrap(~stage, ncol = 2)+
    labs(title = "Biomass allocated to stem")+
    theme(legend.position = "none")
  
  root_plot <- mass_tab%>%
    dplyr::select(week, stage, sp, ends_with("root"))%>%
    ggplot(aes(x = week, y = mean_root, group = factor(sp), colour = factor(sp)))+
    geom_errorbar(aes(min = mean_root - sd_root,
                      max = mean_root + sd_root),
                  stat = "identity")+
    geom_line() +
    geom_point()+
    facet_wrap(~stage, ncol = 2)+
    labs(title = "Biomass allocated to root")+
    theme(legend.position = "none")
  
  repr_plot <- mass_tab%>%
    filter(stage == "a")%>%
    dplyr::select(week, stage, sp, ends_with("repr"))%>%
    ggplot(aes(x = week, y = mean_repr, group = factor(sp), colour = factor(sp)))+
    geom_errorbar(aes(min = mean_repr - sd_repr,
                      max = mean_repr + sd_repr),
                  stat = "identity")+
    geom_line() +
    geom_point()+
    labs(title = "Biomass allocated to reproductive structures")+
    theme(legend.position = "none")
  
  biomass_plot <- plot_grid(leaves_plot, 
                            stem_plot, 
                            root_plot, 
                            repr_plot, 
                            ncol = 2, nrow = 2)
  
  # Total biomass growth curve
  # -------------------------#
  growthcurve_tab <- juvads_allreps_tab%>%
    dplyr::select(week, id, stage, sp, leaves, stem, root, repr)%>%
    filter(stage %in% c("j", "a"))%>%
    group_by(week, stage, sp, id)%>%
    mutate(total = leaves + stem + root + repr)
  
  growthcurve <- growthcurve_tab%>%
    ggplot(aes(x = week, y = total, group = factor(id), colour = factor(id)))+
    geom_line() +
    geom_point()+
    facet_wrap(~stage, ncol = 1,
               scales = "free_y")+
    labs(title = "Individual growth curves")+
    theme(legend.position = "none")
  
  return(list(a = mass_summaryrepli, b = mass_tab, c = growthcurve_tab,
              d = leaves_plot, e = stem_plot, f = repr_plot, g = root_plot,
              h = biomass_plot, i = growthcurve))
}

#' Extract and plot population abundances
popabund <- function(juvads_allreps_tab){
  
  ## extract relevant variables
  pop_tab <- juvads_allreps_tab%>%
    dplyr::select(week,sp,stage,repli)%>%
    ungroup()
  
  ## summarize abundaces/species/week
  spabund_tab <- pop_tab%>%
    group_by(week,sp,repli)%>%
    summarize(abundance = n())%>%
    ungroup()%>%
    group_by(week,sp)%>%
    summarize(mean_abundance = mean(abundance, na.rm = TRUE), 
              sd_abundance = sd(abundance, na.rm = TRUE))%>%
    ungroup()
  
  abund_plot <- ggplot(spabund_tab, 
                          aes(x = week, y = mean_abundance,
			      group = sp, color = factor(sp)))+
    geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
    		      ymax= mean_abundance + sd_abundance), 
                  stat = "identity")+
    geom_line()+
    geom_point(size = 1.25)+
    labs(x = "Year", y = "Abundance (mean +- sd)",
         title = "Species abundance variation")+
    theme(legend.position = "none")
  
  return(list(a = pop_tab, b = spabund_tab, c = abund_plot))
}

#' Extract and plot populations structure
popstruct <- function(pop_tab, seeds_allreps_tab, parentsimID){
  
  ## Population strucuture: absolute abundances
    absltstruct_tab <- pop_tab%>%
      group_by(week, sp, stage)%>%
      summarize(abundance = n())%>%
      ungroup()%>%
      bind_rows(., dplyr::select(seeds_allreps_tab, -mode))%>% #merge offspring file
      ungroup()%>%
      group_by(week, sp, stage)%>%
      summarize(mean_abundance = mean(abundance),
                sd_abundance = sd(abundance))%>%
      ungroup()
    
    ## Population strucuture: relative abundances
    rltvstruct_tab <- pop_tab%>%
      group_by(week, sp, stage)%>%
      summarize(abund = n())%>%
      complete(week, stage, fill = list(abundance = 0))%>%
      group_by(week, sp)%>%
      mutate(total_abund = sum(abund),
             rltv_abund = abund/total_abund)%>%
      complete(week, stage, fill = list(rel_abund = 0))%>%
      ungroup()
  
  ## Graphs
  absltstruct_plot <- ggplot(absltstruct_tab,
                             aes(x = week, y = mean_abundance, color = factor(stage)))+
    geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                      ymax = mean_abundance + sd_abundance),
                  stat = "identity")+
    geom_line()+
    geom_point(size = 1.5)+
    labs(x = "Week", y = "Abundance",
         title = "Population structure")+
    facet_wrap(~sp, ncol = 4, scale = "free_y")+
    theme(legend.position = "bottom")
  
  rltvstruct_plot <- ggplot(rltvstruct_tab,
                            aes(x = week, y = rltv_abund, color = as.factor(stage)))+
    geom_line()+
    geom_point(size = 1.5)+
    labs(x = "Week", y = "Relative abundace",
         title = "Population structure")+
    facet_wrap(~sp, ncol = 4)+
    theme(legend.position = "bottom")

  return(list(a = absltstruct_tab, b = rltvstruct_tab,
              c = absltstruct_plot, d = rltvstruct_plot))
}

#' Calculate mean species richness from replicates and per group of size
spprichness <- function(juvads_allreps_tab, pop_tab, parentsimID, disturbance,tdist){
  
  ## extract richness from output
  spprichness_tab <- pop_tab%>%
    dplyr::select(week, sp, repli)%>%
    group_by(week)%>%
    summarize(richness = length(unique(sp)),
              mean_richness = mean(richness),
              sd_richness = sd(richness))%>%
    ungroup()
  
  ## create plots
  ## identify when disturbance happens, if it does
  spprichness_plot <- ggplot(spprichness_tab,
                             aes(x = week, y = mean_richness))+
    geom_errorbar(aes(ymin = mean_richness - sd_richness,
                      ymax = mean_richness + sd_richness),
                  stat = "identity")+
    geom_line(color = "darkslateblue")+
    geom_point(color = "darkslateblue", size = 1.25)+
    labs(x = "Time", y = "Species richness")
    
    if (disturbance != "n"){
    
        if (disturbance == "a"){
      	   text <- "Area loss"
    	} else if (disturbance == "p"){
      	   text <- "Pollination loss"
    	} else if (disturbance == "ap"){
      	   text <- "Area + pollination loss"
    	}

    spprichness_plot <- spprichness_plot +
                        geom_vline(xintercept = tdist, linetype = 2, color = "red")
      
    }
  
  ## richness per group of size
  groupspprichness_tab  <- juvads_allreps_tab%>%
    dplyr::select(week, sp, seedmass)%>%
    ungroup()%>%
    group_by(week, seedmass)%>%
    summarize(richness = length(unique(sp)))%>%
    ungroup()
  
  groupspprichness_plot <- ggplot(data = groupspprichness_tab, aes(x = week, y = richness))+
    facet_wrap(~seedmass, scales = "free_y", ncol = 1,
               labeller = labeller(seedmass = size_labels))+
    geom_line(color = "darkslateblue")+
    geom_point(color = "darkslateblue", size = 1.25)+
    labs(x = "Time", y = "Species richness")
  
  return(list(a = spprichness_tab, b = spprichness_plot,
              c = groupspprichness_tab, d = groupspprichness_plot))
  
}

#' Calculate species relative abundance and plot species rank-abundance at a given time-step
rankabund <- function(pop_tab, timesteps){
  
  rltvabund_tab <- pop_tab%>%
    filter(week %in% timesteps)%>%
    group_by(week, sp)%>%
    summarize(abundance = n())%>%
    ungroup()%>%
    group_by(sp, week)%>%
    summarize(mean_abundance = mean(abundance))%>%
    ungroup()%>%
    group_by(week)%>%
    mutate(total_abund = sum(mean_abundance),
           rltv_abund = mean_abundance/total_abund)
  
  plot_rankabund <- function(timestep){
    ggplot(filter(rltvabund_tab, week %in% timestep),
           aes(x = reorder(sp, -rltv_abund), y = rltv_abund))+
      geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 80, size = 10, vjust = 0.5))
  }
  
  rankabund_plots <- timesteps%>%
    map(. %>% plot_rankabund)
  
  rankabund_plot <- plot_grid(plotlist = rankabund_plots, 
                              labels = paste("week", timesteps, sep = " "))
  
  return(list(a = rltvabund_tab, b = rankabund_plots, c = rankabund_plot))
}

#' Population structure by group size
#' Plot dynamics of population and total biomass growth per group of size, not species
groupdyn <- function(juvads_allreps_tab, singlestages){
  
  ## create table
  grouppop_tab <- juvads_allreps_tab%>%
    group_by(as.factor(week), sp, seedmass, stage)%>%
    summarize(abundance = n())%>%
    ungroup()%>%
    bind_rows(., dplyr::select(seeds_allreps_tab, -mode))%>% # merge seed info
    group_by(sp)%>% # necessary to fill in seed size according to species
    fill(seedmass)%>%
    ungroup()%>%
    group_by(week, seedmass, stage)%>% 
    summarize(mean_abundance = mean(abundance),
              sd_abundance = sd(abundance))%>%
    ungroup()
  
  groupweight_tab <- juvads_allreps_tab%>%
    dplyr::select(week, id, stage, seedmass, sp, leaves, stem, root, repr)%>%
    mutate(total = leaves + stem + root + repr)%>%
    group_by(week, stage, seedmass)%>%
    summarize(totalmass_mean = mean(total),
              totalmass_sd = sd(total))%>%
    ungroup()
  
  # intra-replicate summary
  groupweight_replitab <- juvads_allreps_tab%>%
    dplyr::select(week, id, stage, seedmass, sp, leaves, stem, root, repr, repli)%>%
    mutate(total = leaves + stem + root + repr)%>%
    group_by(week, stage, seedmass, repli)%>%
    summarize(totalmass_mean = mean(total),
              totalmass_sd = sd(total))%>%
    ungroup()
  
  ## population dynamic
  grouppop_plot <- ggplot(data = grouppop_tab,
                          aes(x = week, y = mean_abundance, colour = stage))+
    geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                      ymax = mean_abundance + sd_abundance),
                  stat = "identity")+
    geom_line()+
    geom_point(size = 1.25)+
    facet_wrap(~seedmass, ncol = 1,
               labeller = labeller(seedmass = size_labels),
               scales = "free_y")+
    labs(title = "Population structure per group size")+
    theme(legend.position = "bottom")
  
  ## plot weigh variation
  groupweight_plot <- ggplot(groupweight_tab,
                             aes(x = week, y = totalmass_mean, colour = stage))+
    geom_errorbar(aes(ymin = totalmass_mean -totalmass_sd, 
                      ymax = totalmass_mean + totalmass_sd))+
    geom_line()+
    geom_point(size = 1.25)+
    facet_wrap(~seedmass, ncol = 1,
               labeller = labeller(seedmass = size_labels),
               scales = "free_y")+
    labs(title = "Total biomass per group size")+
    theme(legend.position = "bottom")
  
  return(list(a = grouppop_tab, b = groupweight_tab, c = groupweight_replitab,
              d = grouppop_plot, e = groupweight_plot))
}

#' Calculate biomass production
production <- function(juvads_allreps_tab){
  
  production_tab <- juvads_allreps_tab%>%
    dplyr::select(week, leaves, stem, root, repr, repli)%>%
    mutate(production = (leaves+stem+root+repr)/(10^3))%>%
    ungroup()%>%
    group_by(week)%>%
    summarize(prod_mean = mean(production),
              prod_sd = sd(production))%>%
    ungroup()
  
  production_replitab <- juvads_allreps_tab%>%
    dplyr::select(week, leaves, stem, root, repr,repli)%>%
    mutate(production = (leaves+stem+root+repr)/(10^3))%>%
    ungroup()%>%
    group_by(week, repli)%>%
    summarize(prod_mean = mean(production))%>%
    ungroup()
  
  production_plot <- ggplot(production_tab, aes(x= week, y = prod_mean))+
    geom_errorbar(aes(ymin = prod_mean - prod_sd,
                      ymax = prod_mean + prod_sd),
                  stat = "identity")+
    geom_point()+
    geom_line()+
    labs(x = "Week",
         y = "Biomass production (kg)")
  
  return(list(a = production_tab, b = production_replitab, c = production_plot))
}

#' Analysis of change in trait values distribution
traitchange  <- function(juvads_allreps_tab, timesteps, species){
  
  if(missing(species)){
    species <- unique(juvads_allreps_tab$sp)
  }
  species <- factor(species)
  
  #Plot trait distribution for the species, at different time steps
  traitvalue_tab <- juvads_allreps_tab%>%
    dplyr::select(-c(kernel, clonality, age, xloc, yloc, leaves, stem, root, repr, mated))
  
  traitvalue4dist <- gather(traitvalue_tab%>%
                              filter(week %in% timesteps),
                            key = trait,
                            value = value,
                            seedmass:bankduration,
                            factor_key = TRUE)%>%
    split(.$trait)
  
  plotdist <- function(spp){
    traitplot <- traitvalue4dist %>%
      map(~ ggplot(.x%>%filter(sp %in% spp),
                   aes(x = sp, y = value))+
            geom_violin()+
            geom_boxplot(width = 0.2)+
            facet_wrap(c("week", "trait"), #facet_grid cannot free y axis
                       nrow = length(unique(timesteps)),
                       scales = "free_y")+
            background_grid(major = "xy", minor = "none"))
    
    return(traitplot)
  }
  
  traitdistribution_plots <- species%>%
    map(. %>% plotdist)
  
  # Plot trait variance over time
  traitsummary_tab  <- juvads_allreps_tab%>%
    dplyr::select(-c(id, stage, kernel, clonality, age, xloc, yloc, 
              leaves, stem, root, repr, mated, repli))%>%
    group_by(week, sp)%>%
    summarize_all(funs(mean = mean,
                       sd = sd))%>%
    ungroup()
  
  
  traitsummary4plot <-  gather(traitsummary_tab,
                               key = trait,
                               value = value,
                               seedmass_mean:bankduration_sd,
                               factor_key = TRUE)%>%
    mutate(trait_metric = str_extract(trait, "^[^_]+"))%>%
    split(.$trait_metric)%>%
    #map(~ dplyr::select(.x, -trait_metric))%>%
    map(~ spread(.x,
                 key = trait,
                 value = value))%>%
    map(~ rename(.x, trait_mean = ends_with("mean"),
                 trait_sd = ends_with("sd")))
  
  plottrait <- function(spp){
    traitplot <- traitsummary4plot %>%
      map(~ ggplot(.x%>%filter(sp %in% spp),
                   aes(x = week, y = trait_mean))+
            geom_line()+
            geom_point(size = 0.5)+
            geom_errorbar(aes(ymin = trait_mean-trait_sd, 
                              ymax = trait_mean+trait_sd), alpha = 0.2)+
            labs(title = as.character(unique(.x$trait_metric))))
    
    return(traitplot)
  }
  
  traitts_plots <- species%>%
    map(. %>% plottrait)    
  
  return(list(a = traitvalue_tab, b = traitdistribution_plots,
              c = traitsummary_tab, d = traitts_plots))
}

#' Analysis of change in trait space
traitspacechange  <- function(traitsdistributions_tab, timesteps){
  
  ## PCA (internal function because it needs to passed to map)
  traitPCA <- function(timestep){
    
    traitpca <- PCA(traitsdistributions_tab%>%
                      dplyr::select(-sp, -stage, -ends_with("_sd"))%>%
                      filter(week %in% timestep)%>%
                      dplyr::select(-week, -id),
                    scale.unit = TRUE,
                    ncp = 5,
                    graph = FALSE)
    
    return(traitpca)
  }
  
  ##separate PCAs for t0 and tend
  traitpcas <- timesteps%>%
    map(., traitPCA)
  
  ## PCA comparing t0 and end
  timepca <- traitPCA(timesteps)
  timepca_plot <- fviz_pca_biplot(timepca,
                                  geom.ind = "point", # show points only (but not "text")
                                  pointshape = 21,
                                  pointsize = 2.5,
                                  fill.ind = factor(dplyr::select(filter(traitsdistributions_tab,
					            week %in% timesteps), week)), # color by time
                                  col.ind = "black",
                                  ##addEllipses = TRUE, # Concentration ellipse,
                                  col.var = factor(c("size", "reprd", "reprd", "reprd", "reprd",
				  	    	     "reprd", "span", "metab", "metab", "metab",
						     "reprd", "size")),
                                  repel = TRUE,
                                  legend.title = list(fill = "Time-step", color = "Traits"))
  
  return(list(a = traitpcas, b = timepca, c = timepca_plot))
  
}

#' Life-history events
#lifehistory <- function(outputsdir){
#  
#  lifeevents_tab <- read_tsv(file.path(outputsdir, "eventslog.txt"), col_names = TRUE)
#  lifeevents_plot <- lifeevents_tab%>%
#    dplyr::select(-age)%>%
#    group_by(week, event, stage)%>%
#    summarize(total = n())%>%
#    ggplot(aes(x = week, y = total, color = event))+
#    geom_line()+
#    facet_wrap(~stage, ncol = 1, scales = "free_y")
#  
#  metabolic_tab <- read_tsv(file.path(outputsdir, "metaboliclog.txt"), col_names = TRUE)
#  metabolic_summary_tab <- metabolic_tab%>%
#    dplyr::select(-age)%>%
#    gather(c(rate, probability), key = "metric", value = value)%>%
#    group_by(stage, event, metric)%>%
#    summarize(mean = mean(value),
#              sd = sd(value))
#  
#  return(list(a = lifeevents_tab, b = lifeevents_plot,
#  	      c = metabolic_tab, d = metabolic_summary_tab))
#}

#' Age distribution of stages
agetraits <- function(juvads_allreps_tab){
  
  meanage_tab <- juvads_allreps_tab%>%
    dplyr::select(-repli)%>%
    group_by(week, sp, stage)%>%
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
    geom_errorbar(aes(ymin=mean_age-sd_age, ymax=mean_age+sd_age))+
    geom_line()
    
  meanage_adtplot <- meanage_tab%>%
    filter(stage == "a")%>%
    ggplot(aes(x=week, y=mean_age, colour=stage))+
    geom_hline(aes(yintercept = mean_span), colour = "grey", alpha = 0.1)+
    geom_vline(aes(xintercept = mean_span), colour = "grey", alpha = 0.1)+
    geom_errorbar(aes(ymin=mean_age-sd_age, ymax=mean_age+sd_age))+
    geom_line()
  
  meanage_plot <- plot_grid(meanage_juvplot, meanage_adtplot,
                            ncol = 2, nrow = 1)
  return(list(a = meanage_tab, b = meanage_plot))
}

##facet_grid(rows = vars(stage), cols = vars(metric), scales = "free_y")+
##scale_color_viridis(discrete = TRUE)

#' Seed production and seed bank
#seeddynamics <- function(juvads_allreps_tab, seeds_allreps_tab){
#  
#  seeddyn_tab <- juvads_allreps_tab%>%
#    filter(stage == "s")%>%
#    dplyr::select(week, sp, stage, repli)%>%
#    group_by(week, sp, stage, repli)%>%
#    summarize(abundance = n())%>%
#    ungroup()%>%
#    mutate(stage = "seedbank")%>%
#    bind_rows(dplyr::select(seeds_allreps_tab, -mode))
#  
#  seeddyn_plot <- ggplot(seeddyn_tab, aes(x = week, y = abundance))+
#    geom_line(aes(group = sp, colour = stage))+
#    geom_point()
#    
#  return(list(a = seeddyn_tab, b = seeddyn_plot))
#}

#### Clean main output ####
#-------------------------#
cleanoutput <- getoutput(parentsimID, repfolder, nreps,
	       	         outputsdir = outputsdir, EDdir = EDdir)  

#### Identify replicates #####
replicates <- orgareplicates(parentsimID, repfolder, nreps)
replicates$a -> juvads_allreps_tab
replicates$b -> seeds_allreps_tab
rm(replicates)

#### Biomasses #####
biomass <- biomass_allocation(juvads_allreps_tab)
biomass$a -> mass_summaryrepli
biomass$b -> mass_tab
biomass$c -> growthcurve_tab
biomass$d -> leaves_plot 
biomass$e -> stem_plot 
biomass$f -> repr_plot
biomass$g -> root_plot
biomass$h -> biomass_plot
biomass$i -> growthcurve_bigplot
rm(biomass)

#### Abundances #####
abund <- popabund(juvads_allreps_tab)
abund$a -> pop_tab
abund$b -> spabund_tab
abund$c -> abund_plot
rm(abund)

#### Population structure #####
population <- popstruct(pop_tab, seeds_allreps_tab, parentsimID)
population$a -> absstruct_tab
population$b -> rltvstruct_tab
population$c -> absstruct_plot
population$d -> rltvstruct_plot
rm(population)

#### Species richness ####
spprich <- spprichness(juvads_allreps_tab, pop_tab, parentsimID, disturbance)
spprich$a -> spprichness_tab
spprich$b -> spprichness_plot 
spprich$c -> groupspprichness_tab
spprich$d -> groupspprichness_plot
rm(spprich)

#### Species rank-abundance ####
## Set up time-steps for which to output derived analysis
#if(!("timesteps" %in% ls())){timesteps <- c(min(pop_tab$week), max(pop_tab$week))}
timesteps <- factor(c(min(pop_tab$week), max(pop_tab$week))) #provided or not, must be converted

rank <- rankabund(pop_tab, timesteps)
rank$a -> rltvabund_tab
rank$b -> rankabund_plots
rank$c -> rankabund_plot
rm(rank)

#### Population structure by group size ####
groups <- groupdyn(juvads_allreps_tab) #specifying singlestages is optional
groups$a -> grouppop_tab
groups$b -> groupweight_tab
groups$c -> groupweight_replitab
groups$d -> grouppop_plot
groups$e -> groupweight_plot
rm(groups)

#### Biomass production ####
prod <- production(juvads_allreps_tab)
prod$a -> production_tab
prod$b -> production_replitab
prod$c -> production_plot
rm(prod)

#### Trait change ####
### trait values
#traitschange <- traitchange(juvads_allreps_tab, timesteps)
#traitschange$a -> traitvalues_tab
#traitschange$b -> traitdistributions_plots
#traitschange$c -> traitssummary_tab
#traitschange$d -> traitts_plots
#rm(traitschange)
### trait space
##traitspace <- traitspacechange(traitsdistributions_tab, timesteps)
##traitspace$a -> traitpcas
##traitspace$b -> timepca
##traitspace$c -> timepca_plot
##rm(traitspace)

#### Life history ####
#lifehistory <- lifehistory(outputsdir)
#lifehistory$a -> events_tab
#lifehistory$b -> events_plot
#lifehistory$c -> metabolic_tab
#lifehistory$d -> metabolic_summary_tab
#rm(lifehistory)

#### Age traits ####
meanage <- agetraits(juvads_allreps_tab)
meanage$a -> meanage_tab
meanage$b -> meanage_plot
#### Seed dynamics ####
#seeddyn <- seeddynamics(juvads_allreps_tab, seeds_allreps_tab)
#seeddyn$a -> seeddyn_tab
#seeddyn$b -> seeddyn_plot
#rm(seeddyn)

#### Save bundle of tabs and plots as RData ####
EDtabs <- objects(name = environment(), all.names = FALSE, pattern = "tab$")
save(list = EDtabs, file = file.path(analysEDdir,
                                     paste(parentsimID, "_tabs.RData", sep = "")))

EDplots <- objects(name = environment(), all.names = FALSE, pattern = "plot$")
save(list = EDplots, file = file.path(analysEDdir,
                                      paste(parentsimID, "_plots.RData", sep = "")))

# Growth curve is rather heavy, dont save an image for it
save(growthcurve_bigplot, file = file.path(analysEDdir,
                                  paste("growthcurve", ".RData", sep = "")))

# Rank-abundance plots are in a list
save(rankabund_plots, file = file.path(analysEDdir,
                                  paste(parentsimID, "rankabunds", ".RData", sep = "")))

# Save lists with plots of trait distribution and values over time
traits_plots <- objects(name = environment(), all.names = FALSE, pattern = "plots$")
save(list = traits_plots, file = file.path(analysEDdir,
                                            paste(parentsimID, "traitsdistributions",
					          ".RData", sep = "")))

# Plot all graphs
map(EDplots, ~ save_plot(filename = file.path(analysEDdir, paste(.x, ".png", sep ="")),
	       	         plot = get(.x)))

#traitvaluesdir <- file.path(analysEDdir, "traitvalues")
#dir.create(traitvaluesdir)

#for(sp in 1:length(unique(juvads_allreps_tab$sp))){
#  for(trait in 1:length(traitdistributions_plots[[sp]])){
#    distributions <- plot_grid(traitdistributions_plots[[sp]][[1]],
#                               traitdistributions_plots[[sp]][[2]],
#                               traitdistributions_plots[[sp]][[3]],
#                               traitdistributions_plots[[sp]][[4]],
#                               traitdistributions_plots[[sp]][[5]],
#                               traitdistributions_plots[[sp]][[6]],
#                               traitdistributions_plots[[sp]][[7]],
#                               traitdistributions_plots[[sp]][[8]],
#                               traitdistributions_plots[[sp]][[9]],
#                               traitdistributions_plots[[sp]][[10]],
#                               nrow = 2,
#                               ncol = 5)
#    ggsave(file = file.path(traitvaluesdir,
#                            paste(unique(juvads_allreps_tab$sp)[sp],"_dist.png", sep = "")),
#           plot = distributions,
#           dpi = 300,
#           height = 30,
#           width = 40,
#           units = "cm")
#    
#    ts <- plot_grid(traitts_plots[[sp]][[1]],
#                   traitts_plots[[sp]][[2]],
#                    traitts_plots[[sp]][[3]],
#                    traitts_plots[[sp]][[4]],
#                    traitts_plots[[sp]][[5]],
#                    traitts_plots[[sp]][[6]],
#                    traitts_plots[[sp]][[7]],
#                    traitts_plots[[sp]][[8]],
#                    traitts_plots[[sp]][[9]],
#                    traitts_plots[[sp]][[10]],
#                    nrow = 2,
#                    ncol = 5)
#    ggsave(file = file.path(traitvaluesdir,
#                            paste(unique(juvads_allreps_tab$sp)[sp],"_ts.png", sep = "")),
#              plot = ts,
#              dpi = 300,
#              height = 30,
#              width = 40,
#              units = "cm")
#  }
#}