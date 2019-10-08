## !Run it from the model's directory!
## This is an R script intended at fast analysis. Details on the RNotebook file of same name
##

## Load packages 
library(tidyverse);
library(viridis);
library(grid);
library(gridExtra);
library(FactoMineR);
library(factoextra); #might require mvtnorm 1.0-6, which does not require R 3.5: devtools::install_version("mvtnorm", version = "1.0-6", repos = "http://cran.us.r-project.org")
library(corrplot);
library(gganimate); #sudo apt-get install cargo install.packages("gifski")
library(cowplot);

## Set directories
source("dirnames_ludmilla.R")

                                        # Set up directory and environement to store analysis
analysEDdir <- file.path(outputsdir, paste(parentsimID, "analysED", sep = "_"))
dir.create(analysEDdir)

## Set theme and colours
theme_set(theme_minimal())
## scale_color_viridis(discrete = TRUE, option = "magma)

## Default values for analysis

#' Organize individual-based output
#' 
#' @param parentsimID The simulation ID
#' @param nreps Number of replicates
#' @param outputsdir The path to the outputs folder, if not default
#' @param EDdir The path to EDM, if not default
getoutput <- function(parentsimID, repfolder, nreps, outputsdir, EDdir = file.path("~/model")){
    
    ## initiliaze list to contain outputs 
    ##outdatalist <- list()
    for(repli in repfolder){
        
        ## get raw outputs
        outraw <- read_tsv(file.path(repli, "orgsweekly.txt"), col_names = TRUE,
                           col_types = cols(week = col_integer(),
                                            id = col_character(),
                                            stage = col_character(),
                                            location = col_character(),
                                            sp = col_character(),
                                            kernel = col_character(),
                                            clonality = col_character(),
                                            seedmass = col_double(),
                                            compartsize = col_double(),
                                            span = col_integer(),
                                            firstflower = col_integer(),
                                            floron = col_integer(),
                                            floroff = col_integer(),
                                            seednumber = col_integer(),
                                            seedon = col_integer(),
                                            seedoff = col_integer(),
                                            bankduration = col_integer(),
                                            b0grow = col_double(),
                                            b0germ = col_double(),
                                            b0mort = col_double(),
                                            age = col_integer(),
                                            leaves = col_double(),
                                            stem = col_double(),
					    root = col_double(),
                                            repr = col_double(),
                                            mated = col_character()));
        ## clean it
        ## take parentheses out of location column ("()")
        loc <- gsub("[\\(|\\)]", "", outraw$location)
        ## split location in 3 columns
        loc <- as.data.frame(matrix(unlist(str_split(loc, ",")), ncol = 2, byrow = T))
        names(loc) = c("xloc", "yloc")
        ## complete and clean table
        outdata <- bind_cols(select(outraw, -location), loc)
        rm(loc)
        ## write single files in its folders
        write.csv(outdata, file.path(repli, paste(parentsimID, "indout.csv", sep = "")), row.names = FALSE)
        ## append to replicates list
        ##outdatalist <- c(outdatalist, 
        ##outdata)
    }
    return(outdata)
}

#' Assemble output files of replicates and identify them
#' Output on juvenile and adult individuals are kept separate because they offspring contains weekly abundances of each species, whereas the juv/adult output is individual-based.
#' @param parentsimID simulation ID
#' @param nreps number of replicates, which identify the folders containing results
orgreplicates <- function(parentsimID, repfolder, nreps){
    
    orgs_complete_tab <- data.frame()
    offspring_complete_tab <- data.frame()
    
    for(i in 1:nreps) {
        sim <- paste(parentsimID, i, sep = "_")
        folder <- repfolder[i]
        ## read outputs of juv/adults and of offspring
        outdatasim <- read.csv(file.path(folder, paste(parentsimID, "indout.csv", sep = "")), header = TRUE)%>%
            mutate(repli = as.factor(rep(sim, nrow(.)))) #include a replication column, to identify replicates
        offspringsim <- read.table(file.path(folder, "offspringproduction.csv"), header = TRUE, sep = "\t")%>%
            mutate(repli = as.factor(rep(sim, nrow(.))))
        
        ## fill in cleanoutput object (necessary step because of bind_rows)
        if(length(orgs_complete_tab) == 0){
            orgs_complete_tab <- outdatasim
            offspring_complete_tab <- offspringsim
        }else{
            ## identify the lines corresponding to it with its replicate ID
            orgs_complete_tab <- bind_rows(orgs_complete_tab, outdatasim, .id = "repli")
            offspring_complete_tab <- bind_rows(offspringsim, offspringsim, .id = "repli")
            offspring_complete_tab$week = factor(offspring_complete_tab$week)
        } 
    }
    return(list(a = orgs_complete_tab, b = offspring_complete_tab))
}

#' Extract and plot species-specific mean and sd of biomass compartments of juveniles and adults.
#' @param outdatalist
#' @param plotit Boolean specifying whether graph shouuld be plotted or not
stagemass <- function(orgs_complete_tab, stages = factor(c("a", "j", "s"))){
    
    ## get biomasses allocated to each compartment
  # summary within replicates
    biomass_tabrepli <- orgs_complete_tab%>%
        select(week, id, stage, sp, leaves, stem, root, repr, repli)%>%
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
    biomass_tab <- orgs_complete_tab%>%
      select(week, id, stage, sp, leaves, stem, root, repr)%>%
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
    leaves_plot <- biomass_tab%>%
      dplyr::select(week, stage, sp, ends_with("leaves"))%>%
      ggplot(aes(x = week, y = mean_leaves, group = factor(sp), colour = factor(sp)))+
      geom_line(position = position_dodge(0.1)) +
      geom_point(position = position_dodge(0.1), size = 0.3)+
      geom_errorbar(aes(min = mean_leaves - sd_leaves,
                        max = mean_leaves + sd_leaves),
                    stat = "identity", colour = "gray50", width = 0.01, position = position_dodge(0.1))+
      facet_wrap(~stage, ncol = 2)+ labs(title = "Biomass allocated to leaves")+
      theme(legend.position = "none")
    
    stem_plot <- biomass_tab%>%
      dplyr::select(week, stage, sp, ends_with("stem"))%>%
      ggplot(aes(x = week, y = mean_stem, group = factor(sp), colour = factor(sp)))+
      geom_line(position = position_dodge(0.1)) +
      geom_point(position = position_dodge(0.1), size = 0.3)+
      geom_errorbar(aes(min = mean_stem - sd_stem,
                        max = mean_stem + sd_stem),
                    stat = "identity", colour = "gray50", width = 0.01, position = position_dodge(0.1))+
      facet_wrap(~stage, ncol = 2)+
      labs(title = "Biomass allocated to stem")+
      theme(legend.position = "none")
    
    root_plot <- biomass_tab%>%
      dplyr::select(week, stage, sp, ends_with("root"))%>%
      ggplot(aes(x = week, y = mean_root, group = factor(sp), colour = factor(sp)))+
      geom_line(position = position_dodge(0.1)) +
      geom_point(position = position_dodge(0.1), size = 0.3)+
      geom_errorbar(aes(min = mean_root - sd_root,
                        max = mean_root + sd_root),
                    stat = "identity", colour = "gray50", width = 0.01, position = position_dodge(0.1))+
      facet_wrap(~stage, ncol = 2)+
      labs(title = "Biomass allocated to root")+
      theme(legend.position = "none") 
    
    repr_plot <- biomass_tab%>%
      filter(stage == "a")%>%
      dplyr::select(week, stage, sp, ends_with("repr"))%>%
      ggplot(aes(x = week, y = mean_repr, group = factor(sp), colour = factor(sp)))+
      geom_line(position = position_dodge(0.1)) +
      geom_point(position = position_dodge(0.1), size = 0.3)+
      geom_errorbar(aes(min = mean_repr - sd_repr,
                        max = mean_repr + sd_repr),
                    stat = "identity", colour = "gray50", width = 0.01, position = position_dodge(0.1))+
      labs(title = "Biomass allocated to reproductive structures")+
      theme(legend.position = "none")
    
    # Total biomass growth curve
    growthcurve_tab <- orgs_complete_tab%>%
      select(week, id, stage, sp, leaves, stem, root, repr)%>%
      filter(stage %in% c("j", "a"))%>%
      group_by(week, stage, sp, id)%>%
      summarize(total = sum(leaves, stem, root, repr))
    
    growthcurve_plot <- growthcurve_tab%>%
      ggplot(aes(x = week, y = total, group = factor(id), colour = factor(id)))+
      geom_line(position = position_dodge(0.1)) +
      geom_point(position = position_dodge(0.1), size = 0.3)+
      facet_wrap(~stage, ncol = 2)+
      scale_colour_viridis(discrete=TRUE)+
      labs(title = "Growth curve of individuals")+
      theme(legend.position = "none")
    
    biomass_plots <- plot_grid(leaves_plot, 
                              stem_plot, 
                              root_plot, 
                              repr_plot, 
                              ncol = 2, nrow = 2)

    return(list(a = vegmass_plottED, b = repmass_plottED, c = biomass_tab))
}

#' Extract and plot population abundances
popabund <- function(orgs_complete_tab){
    
    ## extract relevant variables
    pop_tab <- orgs_complete_tab%>%
        select(week,sp,stage,repli)%>%
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
    
    abund_plottED <- ggplot(spabund_tab, aes(x = week, y = mean_abundance, group = sp, color = factor(sp)))+
        geom_errorbar(aes(ymin = mean_abundance - sd_abundance, ymax= mean_abundance + sd_abundance), 
                      stat = "identity", colour = "gray50", width=.01, position=position_dodge(0.1))+
        geom_line(position=position_dodge(0.1))+
        geom_point(position=position_dodge(0.1))+
        labs(x = "Year", y = "Abundance (mean +- sd)",
             title = "Species abundance variation")+
        scale_color_viridis(discrete = TRUE)+
        theme(legend.position = "none")
    
    return(list(a = pop_tab, b = spabund_tab, c = abund_plottED))
}

#' Extract and plot populations structure
popstruct <- function(pop_tab, offspring_complete_tab, parentsimID, spp){
    
    if (missing(spp)){
        ## Population strucuture: abundances
        weekstruct_tab <- pop_tab%>%
            group_by(week, sp, stage, repli)%>%
            summarize(abundance = n())%>%
            ungroup()%>%
            bind_rows(., select(offspring_complete_tab, -mode))%>% #merge offspring file
            ungroup()%>%
            group_by(week, sp, stage)%>%
            summarize(mean_abundance = mean(abundance),
                      sd_abundance = sd(abundance))
        
        ## Population strucuture: proportions
        relativestruct_tab <- weekstruct_tab%>%
            complete(week, stage, fill = list(abundance = 0))%>%
            group_by(week, sp)%>%
            mutate(proportion = 100*mean_abundance/sum(mean_abundance))%>%
            complete(week, stage, fill = list(proportion = 0))%>%
            ungroup()
    }else{
        spps <- paste("p",spp,sep="-")
        ## Population strucuture: abundances
        weekstruct_tab <- pop_tab%>%
            filter(sp %in% spp)%>%
            group_by(week, sp, stage, repli)%>%
            summarize(abundance = n())%>%
            ungroup()%>%
            bind_rows(., select(offspring_complete_tab, -mode))%>% #merge offspring file
            ungroup()%>%
            group_by(week, sp, stage)%>%
            summarize(mean_abundance = mean(abundance, na.rm = TRUE),
                      sd_abundance = sd(abundance, na.rm = TRUE))
        ## Population strucuture: proportions
        relativestruct_tab <- weekstruct_tab%>%
            filter(sp %in% spps)%>%
            complete(week, stage, fill = list(abundance = 0))%>%
            group_by(week, sp)%>%
            mutate(proportion = mean_abundance/sum(mean_abundance))%>%
            complete(week, stage, fill = list(proportion = 0))%>%
            ungroup()
    }
    ## graphs
    weekstruct_plottED <- ggplot(weekstruct_tab,
                                 aes(x = week, y = mean_abundance, color = factor(stage)))+
        geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                          ymax = mean_abundance + sd_abundance),
                      stat = "identity", colour = "gray50", width=.01, position=position_dodge(0.1))+
        geom_line()+
        geom_point()+
        labs(x = "Week", y = "Abundance",
             title = "Population structure")+
        ##scale_color_discrete("Stages:", labels = c("Adults", "Seeds", "Juveniles"))+
        facet_wrap(~sp, ncol = 4, scale = "free_y")+
        scale_color_viridis(discrete = TRUE)+
        theme(legend.position = "bottom")

	ggsave(filename = file.path(analysEDdir, "weekstruct_plot.png"),
       	       plot = weekstruct_plottED,
               width = 30,
               height = 60,
               unit = "cm",
               dpi = 300)
    
    relativestruct_plottED <- ggplot(relativestruct_tab,
                                     aes(x = week, y = proportion, color = as.factor(stage)))+
        geom_line()+
        geom_point()+
        labs(x = "Week", y = "Relative abundace",
             title = "Population structure")+
        facet_wrap(~sp, ncol = 4)+
        scale_color_viridis(discrete = TRUE)+
        theme(legend.position = "bottom")
    
    ggsave(filename = file.path(analysEDdir, "relativestruct_plot.png"),
       plot = relativestruct_plottED,
       width = 30,
       height = 60,
       unit = "cm",
       dpi = 300)

    return(list(a = weekstruct_tab, b = relativestruct_tab))
}

#' Calculate mean species richness from replicates and per group of size
richness <- function(orgs_complete_tab, pop_tab, parentsimID, disturbance,tdist){
    
    ## extract richness from output
    spprichness_tab <- pop_tab%>%
        group_by(week, repli)%>%
        summarize(richness = length(unique(sp)))%>%
        ungroup()%>%
        group_by(week)%>%
        summarize(mean_richness = mean(richness),
                  sd_richness = sd(richness))%>%
        ungroup()
    
    ## create plots
    ## identify when disturbance happens, if it does
    if(disturbance == "none" && missing(tdist)){
        spprichness_plottED <- ggplot(spprichness_tab,
                                      aes(x = week, y = mean_richness))+
            geom_errorbar(aes(ymin = mean_richness - sd_richness,
                              ymax = mean_richness + sd_richness),
                          stat = "identity", colour = "gray50", width=.01, position=position_dodge(0.1))+
            geom_line(color = "dodgerblue2", size = 1.25)+
            geom_point()+
            scale_color_viridis(discrete = TRUE)
        
    }else{
        
        if (disturbance == "loss"){
            text <- "Area loss"
        } else if (disturbance == "frag"){
            text <- "Fragmentation"
        } else if (disturbance == "poll"){
            text <- "Pollination loss"
        }
        ## TODO: annotate it
        ##my_grob <- grobTree(textGrob(text, x = 0.5, y = 0.9, hjust = 0, 
        ##                             gp = gpar(fontsize = 10, fontface = "italic")))
        spprichness_plottED <- ggplot(spprichness_tab, aes(x = week, y = mean_richness))+
            geom_errorbar(aes(ymin = mean_richness - sd_richness,
                              ymax = mean_richness + sd_richness),
                          stat = "identity", colour = "gray50", width=.01, position=position_dodge(0.1))+
            geom_line(size = 1.25)+
            geom_point()+
            geom_vline(xintercept = tdist, linetype = 2, color = "red")+
            labs(x = "Time", y = "Spp. richness")+
            scale_color_viridis(discrete = TRUE)
        ##annotation_custom(my_grob)+
    }
    
    ## richness per group of size
    groupspprichness_tab  <- orgs_complete_tab%>%
        select(week, sp, seedmass)%>%
        ungroup()%>%
        group_by(week, seedmass)%>%
        summarize(richness = length(unique(sp)))%>%
        ungroup()
    groupspprichness_plottED <- ggplot(data = groupspprichness_tab, aes(x = week, y = richness))+
        facet_wrap(~seedmass, scales = "free_y", ncol = 1)+
        geom_line()
    
    return(list(a = spprichness_tab, b = spprichness_plottED, c = groupspprichness_tab, d = groupspprichness_plottED))
    
}

#' Calculate species relative abundance and plot species rank-abundance at a given time-step
rankabund <- function(pop_tab, timesteps){
    
    relabund_tab <- pop_tab%>%
        filter(week %in% timesteps)%>%
        group_by(week, sp, repli)%>%
        summarize(abundance = n())%>%
        ungroup()%>%
        group_by(sp, week)%>%
        summarize(mean_abundance = mean(abundance, na.rm = TRUE))%>%
        ungroup()%>%
        group_by(week)%>%
        mutate(relabund = mean_abundance/sum(mean_abundance))
    
    plot_rankabund <- function(timestep){
        ggplot(filter(relabund_tab, week %in% timestep),
               aes(x = reorder(sp, -relabund), y = relabund))+
            geom_bar(stat = "identity")+
            scale_color_viridis(discrete = TRUE)+
            theme(axis.text.x = element_text(angle = 50, size = 10, vjust = 0.5))
    }
    
    rankabunds <- timesteps%>%
        map(. %>% plot_rankabund)
    
    return(list(a = relabund_tab, b = rankabunds))
}

#' Population structure by group size
groupdyn <- function(orgs_complete_tab, singlestages){
    
    ## create table
    grouppop_tab <- orgs_complete_tab%>%
        group_by(week, sp, seedmass, stage, repli)%>%
        summarize(abundance = n())%>%
        ungroup()%>%
        bind_rows(., select(offspring_complete_tab, -mode))%>% # merge seed info
        group_by(sp)%>% # necessary to fill in seed size according to species
        fill(seedmass)%>%
        ungroup()%>%
        group_by(week, seedmass, stage)%>% 
        summarize(mean_abundance = mean(abundance),
                  sd_abundance = sd(abundance))%>%
        ungroup()
    
    ## create plots 
    ## filter data, if necessary
    if (missing(singlestages)){
        popdata_to_plot <- grouppop_tab
        weightdata_to_plot <- orgs_complete_tab
    }else{
        popdata_to_plot <- grouppop_tab %>%
            filter(stage %in% singlestages)
        weightdata_to_plot <- orgs_complete_tab%>%
            filter(stage %in% singlestages)
    }
    ## plot it
    grouppop_plot <- ggplot(data = popdata_to_plot,
                            aes(x = week, y = mean_abundance, colour = stage))+
        geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                          ymax = mean_abundance + sd_abundance),
                      stat = "identity", position = position_dodge(0.1), colour = "gray50", width = 0.01)+
        geom_line()+
        geom_point()+
        facet_wrap(~seedmass, ncol = 1, nrow = 3)+
        labs(title = "Population structure per group size")+
        scale_color_viridis(discrete = TRUE)
    
    ## plot weigh variation
    groupweight_plot <- ggplot(data = weightdata_to_plot%>%
                                   group_by(week, seedmass, stage, repli)%>%
                                   summarize(repli_mean_weight = mean(veg, na.rm = TRUE))%>%
                                   ungroup()%>%
                                   group_by(week, seedmass, stage)%>%
                                   summarize(mean_weight = mean(repli_mean_weight, na.rm = TRUE), 
                                             sd_weight = sd(repli_mean_weight, na.rm = TRUE)),
                               aes(x = week, y = mean_weight, colour = stage))+
        geom_errorbar(aes(ymin = mean_weight-sd_weight, ymax = mean_weight + sd_weight), 
                      colour = "black", width=.01, position = position_dodge(0.1))+
        geom_line(position = position_dodge(0.1))+
        geom_point(position = position_dodge(0.1))+
        scale_color_viridis(discrete = TRUE)+
        facet_wrap(~seedmass, ncol = 1, nrow = 3)
    
    return(list(a = grouppop_tab, b = grouppop_plot, c = groupweight_plot))
    
}

#' Calculate biomass production
production <- function(orgs_complete_tab){
    
    production_tab <- orgs_complete_tab%>%
        select(week,veg,repr,repli)%>%
        group_by(week, repli)%>%
        summarize(production = sum(veg,repr)/(10^3))%>%
        ungroup()%>%
        group_by(week)%>%
        summarize(mean_prod = mean(production, na.rm = TRUE),
                  sd_prod = sd(production, na.rm = TRUE))%>%
        ungroup()
    
    production_plot<- ggplot(production_tab, aes(x= week, y = mean_prod))+
        geom_errorbar(aes(ymin = mean_prod - sd_prod,
                          ymax = mean_prod + sd_prod),
                      stat = "identity", position = position_dodge(0.01), colour = "gray50", width = 0.01)+
        geom_point()+
        geom_line()+
        labs(x = "Week",
             y = "Biomass production (kg)")+
        scale_color_viridis(discrete = TRUE)
    
    return(production_plot)
}

#' Analysis of change in trait values distribution

traitchange  <- function(orgs_complete_tab, timesteps, species){
    
    if(missing(species)){
        species <- unique(orgs_complete_tab$sp)
    }
    species <- factor(species)
    
                                        #Plot trait distribution for the species, at different time steps
    traitvalue_tab <- orgs_complete_tab%>%
        select(-c(kernel, clonality, age, xloc, yloc, veg, repr, mated))

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
                    ##geom_dotplot(color = "grey31", fill = "white", alpha = 0.8)+
                    geom_boxplot(width = 0.2)+
                    facet_wrap(c("week", "trait"), #facet_grid cannot free y axis
                               nrow = length(unique(timesteps)),
                               scales = "free_y")+
                    background_grid(major = "xy", minor = "none")+
                    scale_color_viridis(discrete = TRUE))

        return(traitplot)
    }
    
    traitdistribution_plots <- species%>%
        map(. %>% plotdist)

                                        # Plot trait variance over time
    traitsummary_tab  <- orgs_complete_tab%>%
        select(-c(id, stage, kernel, clonality, age, xloc, yloc, veg, repr, mated, b0mort, b0germ, b0grow, repli))%>%
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
                                        #map(~ select(.x, -trait_metric))%>%
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
                                      ymax = trait_mean+trait_sd),
                                  colour = "gray", alpha = 0.2)+
                    scale_color_viridis(discrete = TRUE)+
                    labs(title = as.character(unique(.x$trait_metric))))
        
        return(traitplot)
    }

    traitts_plots <- species%>%
        map(. %>% plottrait)    

    return(list(a = traitvalue_tab, b = traitdistribution_plots, c = traitsummary_tab, d = traitts_plots))
}

#' Analysis of change in trait space

traitspacechange  <- function(traitsdistributions_tab, timesteps){
    
    ## PCA (internal function because it needs to passed to map)
    traitPCA <- function(timestep){
        
        traitpca <- PCA(traitsdistributions_tab%>%
                        select(-sp, -stage, -ends_with("_sd"))%>%
                        filter(week %in% timestep)%>%
                        select(-week, -id),
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
                                    fill.ind = factor(select(filter(traitsdistributions_tab, week %in% timesteps), week)), # color by time
                                    col.ind = "black",
                                    ##addEllipses = TRUE, # Concentration ellipse,
                                    col.var = factor(c("size", "reprd", "reprd", "reprd", "reprd", "reprd", "span", "metab", "metab", "metab","reprd", "size")),
                                    repel = TRUE,
                                    legend.title = list(fill = "Time-step", color = "Traits"))
    
    return(list(a = traitpcas, b = timepca, c = timepca_plot))
    
}

#' Life-history events
lifehistory <- function(outputsdir){
    
    lifeevents_tab <- read_tsv(file.path(outputsdir, "eventslog.txt"), col_names = TRUE)
    lifeevents_plot <- lifeevents_tab%>%
        select(-age)%>%
        group_by(week, event, stage)%>%
        summarize(total = n())%>%
        ggplot(aes(x = week, y = total, color = event))+
        geom_line()+
        facet_wrap(~stage, ncol = 1, scales = "free_y")+
        scale_color_viridis(discrete = TRUE)

    metabolic_tab <- read_tsv(file.path(outputsdir, "metaboliclog.txt"), col_names = TRUE)
    metabolic_summary_tab <- metabolic_tab%>%
        select(-age)%>%
        gather(c(rate, probability), key = "metric", value = value)%>%
        group_by(stage, event, metric)%>%
        summarize(mean = mean(value),
                  sd = sd(value))

    return(list(a = lifeevents_tab, b = lifeevents_plot, c = metabolic_tab, d = metabolic_summary_tab))
}

#' Age distribution of stages
agetraits <- function(orgs_complete_tab){

    meanage_plot <- orgs_complete_tab%>%
        select(-repli)%>%
        group_by(week, sp, stage)%>%
        summarize(mean_age = mean(age),
                  sd_age = sd(age),
                  mean_span = mean(span),
                  sd_span = sd(span),
                  mean_firstflower = mean(firstflower),
                  sd_firstflower = sd(firstflower))%>%
        ungroup()%>%
        ggplot(aes(x=week, y=mean_age, colour=stage))+
        geom_line()+
        geom_errorbar(aes(ymin=mean_age-sd_age, ymax=mean_age+sd_age))+
        geom_hline(aes(yintercept = mean_span), colour = "grey", alpha = 0.1)+
        geom_hline(aes(yintercept = mean_firstflower), colour = "gold", alpha = 0.1)+
        geom_vline(aes(xintercept = mean_span), colour = "grey", alpha = 0.1)+
        geom_vline(aes(xintercept = mean_firstflower), colour = "gold", alpha = 0.1)+
        scale_color_viridis(discrete = TRUE)+
        theme(legend.position = "bottom")

    return(meanage_plot)
}

##facet_grid(rows = vars(stage), cols = vars(metric), scales = "free_y")+
##scale_color_viridis(discrete = TRUE)

#' Seed production and seed bank
seeddynamics <- function(orgs_complete_tab, offspring_complete_tab){

seeddyn_tab <- orgs_complete_tab%>%
  filter(stage == "s")%>%
  select(week, sp, stage, repli)%>%
  group_by(week, sp, stage, repli)%>%
  summarize(abundance = n())%>%
  ungroup()%>%
  mutate(stage = "seedbank")%>%
  bind_rows(select(offspring_complete_tab, -mode))
  
seeddyn_plot <- ggplot(seeddyn_tab, aes(x = week, y = abundance))+
  geom_point(size = 0.5)+
  geom_line(aes(group = sp, colour = stage))+
  scale_colour_viridis(discrete = TRUE)+
  theme(legend.position = "bottom")

  return(list(a = seeddyn_tab, b = seeddyn_plot))
}
############################################################################
##                         Organize analysis output                       ##
############################################################################

cleanoutput <- getoutput(parentsimID, repfolder, nreps, outputsdir = outputsdir, EDdir = EDdir)  

## Identify replicates
replicates <- orgreplicates(parentsimID, repfolder, nreps)
replicates$a -> orgs_complete_tab
replicates$b -> offspring_complete_tab
rm(replicates)

## Individual vegetative and reproductive biomasses of juveniles and adults (NOT seeds)
#mass <- stagemass(orgs_complete_tab,TRUE) # spp is an optional argument and stage is set to default
#mass$a -> vegmass_plot 
#mass$b -> repmass_plot 
#mass$c -> biomass_tab
#rm(mass)

## Species abundance variation
abund <- popabund(orgs_complete_tab)
abund$a -> pop_tab
abund$b -> spabund_tab
abund$c -> abund_plot
rm(abund)

## Population structure variation
population <- popstruct(pop_tab, offspring_complete_tab, parentsimID) #specifying species is optional
population$a -> weekstruct_tab
population$b -> weekstruct_plot
population$c -> relativestruct_tab
population$d -> relativestruct_plot
rm(population)

## Species richness
##rich <- richness(pop_tab, parentsimID, disturbance) # tdist is optional
##rich$a -> spprichness_tab 
##rich$b -> spprichness_plot
#rich$c -> groupspprichness_tab 
#rich$d -> groupspprichness_plot
#rm(rich)

## Set up time-steps for which to output derived analysis
                                        #if(!("timesteps" %in% ls())){timesteps <- c(min(pop_tab$week), max(pop_tab$week))}
timesteps <- factor(c(min(pop_tab$week), max(pop_tab$week))) #provided or not, must be converted

## Species rank-abundance 
rank <- rankabund(pop_tab, timesteps)
rank$a -> relabund_tab
rank$b -> rankabunds
rm(rank)

## Population structure by group size
groups <- groupdyn(orgs_complete_tab) #specifying singlestages is optional
groups$a -> grouppop_tab
groups$b -> grouppop_plot
groups$c -> groupweight_plot
rm(groups)

## Biomass production
production_plot <- production(orgs_complete_tab)

## Trait change
### trait values
traitschange <- traitchange(orgs_complete_tab, timesteps)
traitschange$a -> traitvalues_tab
traitschange$b -> traitdistributions_plots
traitschange$c -> traitssummary_tab
traitschange$d -> traitts_plots
rm(traitschange)
### trait space
##traitspace <- traitspacechange(traitsdistributions_tab, timesteps)
##traitspace$a -> traitpcas
##traitspace$b -> timepca
##traitspace$c -> timepca_plot
##rm(traitspace)

## Life history
lifehistory <- lifehistory(outputsdir)
lifehistory$a -> events_tab
lifehistory$b -> events_plot
lifehistory$c -> metabolic_tab
lifehistory$d -> metabolic_summary_tab

## Age traits
meanage_plot <- agetraits(orgs_complete_tab)

## Seed dynamics
seeddyn <- seeddynamics(orgs_complete_tab, offspring_complete_tab)
seeddyn$a -> seeddyn_tab
seeddyn$b -> seeddyn_plot

## Save bundle of tabs and plots as RData
EDtabs <- objects(name = environment(), all.names = FALSE, pattern = "_tab$")
save(list = EDtabs, file = file.path(analysEDdir,
                                     paste(parentsimID, "_tabs.RData", sep = "")))

EDplots <- objects(name = environment(), all.names = FALSE, pattern = "_plot$")
save(list = EDplots, file = file.path(analysEDdir,
                                      paste(parentsimID, "_plots.RData", sep = "")))

## Rank-abundance plots are in a list
save(rankabunds, file = file.path(analysEDdir,
                                                     paste(parentsimID, "rankabunds", ".RData", sep = "")))

## Save lists with plots of trait distribution and values over time
traitplots <- objects(name = environment(), all.names = FALSE, pattern = "_plots$")
save(list = traitplots, file = file.path(analysEDdir,
                                         paste(parentsimID, "traitsdistributions", ".RData", sep = "")))

## Plot all graphs
map(EDplots, ~ save_plot(filename = file.path(analysEDdir, paste(.x, ".png", sep ="")), plot = get(.x)))

traitvaluesdir <- file.path(analysEDdir, "traitvalues")
dir.create(traitvaluesdir)

for(sp in 1:length(unique(orgs_complete_tab$sp))){
  for(trait in 1:length(traitdistributions_plots[[sp]])){
    distributions <- plot_grid(traitdistributions_plots[[sp]][[1]],
              traitdistributions_plots[[sp]][[2]],
              traitdistributions_plots[[sp]][[3]],
              traitdistributions_plots[[sp]][[4]],
              traitdistributions_plots[[sp]][[5]],
              traitdistributions_plots[[sp]][[6]],
              traitdistributions_plots[[sp]][[7]],
              traitdistributions_plots[[sp]][[8]],
              traitdistributions_plots[[sp]][[9]],
              traitdistributions_plots[[sp]][[10]],
              nrow = 2,
              ncol = 5)
      ggsave(file = file.path(traitvaluesdir,
                              paste(unique(orgs_complete_tab$sp)[sp],"_dist.png", sep = "")),
             plot = distributions,
             device = "png",
             dpi = 300)
    ts <- plot_grid(traitts_plots[[sp]][[1]],
              traitts_plots[[sp]][[2]],
              traitts_plots[[sp]][[3]],
              traitts_plots[[sp]][[4]],
              traitts_plots[[sp]][[5]],
              traitts_plots[[sp]][[6]],
              traitts_plots[[sp]][[7]],
              traitts_plots[[sp]][[8]],
              traitts_plots[[sp]][[9]],
              traitts_plots[[sp]][[10]],
              nrow = 2,
              ncol = 5)
      ggsave(file = file.path(traitvaluesdir,
                              paste(unique(orgs_complete_tab$sp)[sp],"_ts.png", sep = "")),
             plot = ts,
             device = "png",
             dpi = 300)
  }
}