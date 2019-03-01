# !Run it from the model's directory!
# This is an R script intended at faster analysis. Details on the RNotebook file of same name

# Set up directory and environement to store analysis
analysEDdir <- file.path(outdir, parentsimID, paste(parentsimID, "analysED", sep = "_"))
dir.create(analysEDdir)

# Load packages
library(tidyverse);
library(viridis);
library(grid);
library(gridExtra);

#' Organize individual-based output
#' 
#' @param parentsimID The simulation ID
#' @param nreps Number of replicates
#' @param outdir The path to the outputs folder, if not default
#' @param EDdir The path to EDM, if not default
getoutput <- function(parentsimID, nreps, outdir = file.path("~/model/EDoutputs"), EDdir = file.path("~/model")){
  
  # initiliaze list to contain outputs 
  #outdatalist <- list()
  for(rep in nreps){
    repfolder <- file.path(outdir, parentsimID, paste(parentsimID, rep, sep = "_")); #parentsimID
    
    # get raw outputs
    outraw <- as_tibble(read.table(file.path(repfolder, "orgsweekly.txt"), header = TRUE, sep = "\t"));
    # clean it
    ## take parentheses out of location column ("()")
    loc <- gsub("[\\(|\\)]", "", outraw$location)
    ## split location in 3 columns
    loc <- as.data.frame(matrix(unlist(str_split(loc, ",")), ncol = 2, byrow = T))
    names(loc) = c("xloc", "yloc")
    ## complete and clean table
    outdata <- as_tibble(cbind(select(outraw, -one_of("location")), loc))
    rm(loc)
    # write single files in its folders
    write.csv(outdata, file.path(repfolder, paste(parentsimID, "indout.csv", sep = "")), row.names = FALSE)
    # append to replicates list
    #outdatalist <- c(outdatalist, 
    #outdata)
  }
  return(outdata)
}

#' Assemble output files of replicates and identify them
#' Output on juvenile and adult individuals are kept separate because they offspring contains weekly abundances of each species, whereas the juv/adult output is individual-based.
#' @param parentsimID simulation ID
#' @param nreps number of replicates, which identify the folders containing results
organizereplicates <- function(parentsimID, nreps){
  
  output_replicates <- numeric(0)
  offspring_replicates <- numeric(0)
  
  for(sim in c(paste(parentsimID, seq(1, nreps), sep = "_"))) {
    repfolder <- file.path(outdir, parentsimID, sim); # replicate folder
    # read outputs of juv/adults and of offspring
    outdatasim <- read.csv(file.path(repfolder, paste(parentsimID, "indout.csv", sep = "")), header = TRUE)%>%
      mutate(repli = as.factor(rep(sim, nrow(.)))) #include a replication column, to identify replicates
    offspringsim <- read.table(file.path(repfolder, "offspringproduction.csv"), header = TRUE, sep = "\t")%>%
      mutate(repli = as.factor(rep(sim, nrow(.))))
    
    # fill in cleanoutput object (necessary step because of bind_rows)
    if(length(output_replicates) == 0){
      output_replicates <- outdatasim
      offspring_replicates <- offspringsim
    }else{
      # identify the lines corresponding to it with its replicate ID
      output_replicates <- bind_rows(output_replicates, outdatasim, .id = "repli")
      offspring_replicates <- bind_rows(offspringsim, offspringsim, .id = "repli")
    } 
  }
  return(list(a = output_replicates, b = offspring_replicates))
}

#' Extract and plot species-specific mean and sd biomass compartiments (these are means of each simulation' means   
#' @param outdatalist
#' @param plotit Boolean specifying whether graph shouuld be plotted or not
stagemass <- function(output_replicates, stages = factor(c("a", "j", "s"))){
  
  # get biomasses
  biomass_tab <- output_replicates%>%
    select(week, id, stage, sp, veg, repr, repli)%>%
    group_by(week, stage, sp, repli)%>%
    summarize(repli_mean_reprmass = mean(repr),
              repli_mean_vegmass = mean(veg))%>%
    ungroup()%>%
    group_by(week, stage, sp)%>%
    summarize(mean_reprmass = mean(repli_mean_reprmass),
              mean_vegmass = mean(repli_mean_vegmass),
              sd_reprmass = sd(repli_mean_reprmass),
              sd_vegmass = sd(repli_mean_vegmass))%>%
    ungroup()
  
  # Filter all adult individuals and use colours to distinguish the lines representing each individual
  vegmass_plottED <- ggplot(filter(biomass_tab, stage %in% stages),
                            aes(x=week, y= mean_vegmass, group = factor(stage), color = factor(stage)))+
    geom_line(position = position_dodge(0.1)) +
    geom_point(position = position_dodge(0.1))+
    geom_errorbar(aes(min = mean_vegmass - sd_vegmass,
                      max = mean_vegmass + sd_vegmass),
                  stat = "identity", colour = "gray50", width = 0.01, position = position_dodge(0.1))+
    facet_grid(factor(sp) ~ factor(stage))+
    labs(x = "Week", y = "Compartiment biomass (g)",
         title = "Weekly vegetative biomass")+
    theme(legend.position = "none")
  
  repmass_plottED <- ggplot(filter(biomass_tab, stage %in% stages),
                            aes(x=week, y= mean_reprmass, group = factor(stage), color = factor(stage)))+
    geom_line(position = position_dodge(0.1)) +
    geom_point(position = position_dodge(0.1))+
    geom_errorbar(aes(min = mean_reprmass - sd_reprmass,
                      max = mean_reprmass + sd_reprmass),
                  stat = "identity", colour = "gray50", width = 0.01, position = position_dodge(0.1))+
    labs(x = "Week", y = "Compartiment biomass (g)",
         title = "Weekly reproductive biomass")+
    facet_grid(factor(sp) ~ factor(stage))+
    theme(legend.position = "none")
  
  return(list(a = vegmass_plottED, b = repmass_plottED, c = biomass_tab))
}

#' Extract and plot population abundances
popabund <- function(output_replicates){
  
  # extract relevant variables
  pop_tab <- output_replicates%>%
    select(week,sp,stage,repli)%>%
    ungroup()
  
  # summarize abundaces/species/week
  spabund_tab <- pop_tab%>%
    group_by(week,sp,repli)%>%
    summarize(abundance = n())%>%
    ungroup()%>%
    group_by(week,sp)%>%
    summarize(mean_abundance = mean(abundance), sd_abundance = sd(abundance))%>%
    ungroup()
  
  abund_plottED <- ggplot(spabund_tab, aes(x = week, y = mean_abundance, group = sp, color = factor(sp)))+
    geom_errorbar(aes(ymin = mean_abundance - sd_abundance, ymax= mean_abundance + sd_abundance), 
                  stat = "identity", colour = "gray50", width=.01, position=position_dodge(0.1))+
    geom_line(position=position_dodge(0.1))+
    geom_point(position=position_dodge(0.1))+
    labs(x = "Year", y = "Abundance (mean +- sd)",
         title = "Species abundance variation")+
    theme(legend.position = "none")+
    theme_minimal()
  
  return(list(a = pop_tab, b = spabund_tab, c = abund_plottED))
}

#' Extract and plot populations structure
popstruct <- function(pop_tab, offspring_replicates, parentsimID, spp){
  
  if (missing(spp)){
    # Population strucuture: abundances
    weekstruct_tab <- pop_tab%>%
      group_by(week, sp, stage, repli)%>%
      summarize(abundance = n())%>%
      bind_rows(., select(offspring_replicates, -mode))%>% #merge offspring file
      ungroup()%>%
      group_by(week, sp, stage)%>%
      summarize(mean_abundance = mean(abundance),
                sd_abundance = sd(abundance))
    
    # Population strucuture: proportions
    relativestruct_tab <- weekstruct_tab%>%
      complete(week, stage, fill = list(abundance = 0))%>%
      group_by(week, sp)%>%
      mutate(proportion = 100*mean_abundance/sum(mean_abundance))%>%
      complete(week, stage, fill = list(proportion = 0))%>%
      ungroup()
  }else{
    spps <- paste("p",spp,sep="-")
    # Population strucuture: abundances
    weekstruct_tab <- pop_tab%>%
      filter(sp %in% spp)%>%
      group_by(week, sp, stage, repli)%>%
      summarize(abundance = n())%>%
      bind_rows(., select(offspring_replicates, -mode))%>% #merge offspring file
      ungroup()%>%
      group_by(week, sp, stage)%>%
      summarize(mean_abundance = mean(abundance),
                sd_abundance = sd(abundance))
    # Population strucuture: proportions
    relativestruct_tab <- weekstruct_tab%>%
      filter(sp %in% spps)%>%
      complete(week, stage, fill = list(abundance = 0))%>%
      group_by(week, sp)%>%
      mutate(proportion = mean_abundance/sum(mean_abundance))%>%
      complete(week, stage, fill = list(proportion = 0))%>%
      ungroup()
  }
  # graphs
  weekstruct_plottED <- ggplot(weekstruct_tab,
                               aes(x = week, y = mean_abundance, color = factor(stage)))+
    geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                      ymax = mean_abundance + sd_abundance),
                  stat = "identity", colour = "gray50", width=.01, position=position_dodge(0.1))+
    geom_line()+
    geom_point()+
    labs(x = "Week", y = "Abundance",
         title = "Population structure")+
    #scale_color_discrete("Stages:", labels = c("Adults", "Seeds", "Juveniles"))+
    facet_wrap(~sp, nrow = length(unique(weekstruct_tab$sp)))+
    theme(legend.position = "none")+
    theme_minimal()
  
  relativestruct_plottED <- ggplot(relativestruct_tab,
                                   aes(x = week, y = proportion, color = as.factor(stage)))+
    geom_line()+
    geom_point()+
    labs(x = "Week", y = "Relative abundace",
         title = "Population structure")+
    facet_wrap(~sp, nrow = length(unique(weekstruct_tab$sp)))+
    theme(legend.position = "none")+
    theme_minimal()
  
  return(list(a = weekstruct_tab, b = weekstruct_plottED,
              c = relativestruct_tab, d = relativestruct_plottED))
}

#' Calculate mean species richness from replicates
richness <- function(pop_tab,parentsimID,disturbance,tdist){
  
  # extract richness from output
  spprichness_tab <- pop_tab%>%
    group_by(week, repli)%>%
    summarize(richness = length(unique(sp)))%>%
    ungroup()%>%
    group_by(week)%>%
    summarize(mean_richness = mean(richness),
              sd_richness = sd(richness))%>%
    ungroup()
  
  # create plots
  ## identify when disturbance happens, if it does
  if (disturbance == "none" && missing(tdist)){
    spprichness_plottED <- ggplot(spprichness_tab,
                                  aes(x = week, y = mean_richness))+
      geom_errorbar(aes(ymin = mean_richness - sd_richness,
                        ymax = mean_richness + sd_richness),
                    stat = "identity", colour = "gray50", width=.01, position=position_dodge(0.1))+
      geom_line(color = "dodgerblue2", size = 1.25)+
      geom_point()+
      theme_minimal()
  } else {
    
    if (disturbance == "loss"){
      text <- "Area loss"
    } else if (disturbance == "frag"){
      text <- "Fragmentation"
    } else if (disturbance == "poll"){
      text <- "Pollination loss"
    }
    # TODO: annotate it
    #my_grob <- grobTree(textGrob(text, x = 0.5, y = 0.9, hjust = 0, 
    #                             gp = gpar(fontsize = 10, fontface = "italic")))
    spprichness_plottED <- ggplot(spprichness_tab, aes(x = week, y = mean_richness))+
      geom_errorbar(aes(yim = mean_richness - sd_richness,
                        ymax = mean_richness + sd_richness),
                    stat = "identity", colour = "gray50", width=.01, position=position_dodge(0.1))+
      geom_line(color = "dodgerblue2", size = 1.25)+
      geom_point()+
      geom_vline(xintercept = tdist, linetype = 2, color = "red")+
      labs(x = "Time", y = "Spp. richness")+
      theme_minimal()
    #annotation_custom(my_grob)+
  }
  
  return(list(a = spprichness_tab, b = spprichness_plottED))
  
}

#' Calculate species relative abundance and plot species rank-abundance at a given time-step
rankabund <- function(pop_tab, timestep){
  
  if(!(timestep %in% pop_tab$week)){stop("\'timestep\' for calculating the rank abundance is not available in the outpout")}
  # get relative abundances
  relabund_tab <- pop_tab%>%
    filter(week == timestep)%>%
    group_by(sp, repli)%>%
    summarize(abundance = n())%>%
    ungroup()%>%
    group_by(sp)%>%
    summarize(mean_abundance = mean(abundance))%>%
    ungroup()%>%
    mutate(relabund = mean_abundance/sum(mean_abundance))
  # plot it    
  rankabund_plottED <- ggplot(relabund_tab,
                              aes(x = reorder(sp, -relabund), y = relabund))+
    geom_bar(stat = "identity")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 50, size = 10, vjust = 0.5))
  
  return(list(a = relabund_tab, b = rankabund_plottED))
}

#' Population structure by group size
groupdyn <- function(output_replicates, singlestages){
  
  # create table
  grouppop_tab <- output_replicates%>%
    group_by(week, sp, e_mu, stage, repli)%>%
    summarize(abundance = n())%>%
    ungroup()%>%
    bind_rows(., select(offspring_replicates, -mode))%>% # merge seed info
    group_by(sp)%>% # necessary to fill in seed size according to species
    fill(e_mu)%>%
    ungroup()%>%
    group_by(week, e_mu, stage)%>% 
    summarize(mean_abundance = mean(abundance),
              sd_abundance = sd(abundance))%>%
    ungroup()
  
  # create plots 
  ## filter data, if necessary
  if (missing(singlestages)){
    popdata_to_plot <- grouppop_tab
    weightdata_to_plot <- output_replicates
  }else{
    popdata_to_plot <- grouppop_tab %>%
      filter(stage %in% singlestages)
    weightdata_to_plot <- output_replicates%>%
      filter(stage %in% singlestages)
  }
  ## plot it
  grouppop_plottED <- ggplot(data = popdata_to_plot,
                             aes(x = week, y = mean_abundance, colour = stage))+
    geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                      ymax = mean_abundance + sd_abundance),
                  stat = "identity", position = position_dodge(0.1), colour = "gray50", width = 0.01)+
    geom_line()+
    geom_point()+
    facet_wrap(~e_mu, ncol = 1, nrow = 3)+
    labs(title = "Population structure per group size")+
    theme_minimal()
  
  # plot weigh variation
  groupweight_plottED <- ggplot(data = weightdata_to_plot%>%
                                  group_by(week, e_mu, stage, repli)%>%
                                  summarize(repli_mean_weight = mean(veg))%>%
                                  ungroup()%>%
                                  group_by(week, e_mu, stage)%>%
                                  summarize(mean_weight = mean(repli_mean_weight), 
                                            sd_weight = sd(repli_mean_weight)),
                                aes(x = week, y = mean_weight, colour = stage))+
    geom_errorbar(aes(ymin = mean_weight-sd_weight, ymax = mean_weight + sd_weight), 
                  colour = "black", width=.01, position = position_dodge(0.1))+
    geom_line(position = position_dodge(0.1))+
    geom_point(position = position_dodge(0.1))+
    facet_wrap(~e_mu, ncol = 1, nrow = 3)+
    theme_bw()
  
  return(list(a = grouppop_tab, b = grouppop_plottED, c = groupweight_plottED))
  
}

#' Calculate biomass production
production <- function(output_replicates){
  
  production_tab <- output_replicates%>%
    select(week,veg,repr,repli)%>%
    group_by(week, repli)%>%
    summarize(production = sum(veg,repr)/(10^3))%>%
    ungroup()%>%
    group_by(week)%>%
    summarize(mean_prod = mean(production),
              sd_prod = sd(production))%>%
    ungroup()
  
  production_plottED<- ggplot(production_tab, aes(x= week, y = mean_prod))+
    geom_errorbar(aes(ymin = mean_prod - sd_prod,
                      ymax = mean_prod + sd_prod),
                  stat = "identity", position = position_dodge(0.01), colour = "gray50", width = 0.01)+
    geom_point()+
    geom_line()+
    labs(x = "Week",
         y = "Biomass production (kg)")+
    theme_minimal()
  
  return(production_plottED)
}

# Analysis
cleanoutput <- getoutput(parentsimID, nreps, outdir = outdir, EDdir = EDdir)  

## Identify replicates
replicates <- organizereplicates(parentsimID,nreps)
replicates$a -> output_replicates
replicates$b -> offspring_replicates
rm(replicates)

## Individual vegetative and reproductive biomasses of juveniles and adults (NOT seeds)
mass <- stagemass(output_replicates,TRUE) # spp is an optional argument and stage is set to default
mass$a -> vegmass_plottED 
mass$b -> repmass_plottED 
mass$c -> biomass_tab
rm(mass)

## Species abundance variation
abund <- popabund(output_replicates)
abund$a -> pop_tab
abund$b -> spabund_tab
abund$c -> abund_plottED
rm(abund)

## Population structure variation
pop <- popstruct(pop_tab, offspring_replicates, parentsimID) #specifying species is optional
pop$a -> weekstruct_tab
pop$b -> weekstruct_plottED
pop$c -> relativestruct_tab
pop$d -> relativestruct_plottED
rm(pop)

## Species richness
rich <- richness(pop_tab, parentsimID, disturbance) # tdist is optional
rich$a -> spprichness_tab 
rich$b -> spprichness_plottED
rm(rich)

## Species rank-abundance 
rank <- rankabund(pop_tab, timestep = 13)
rank$a -> relabund_tab
rank$b -> rankabund_plottED
rm(rank)

## Population structure by group size
groups <- groupdyn(output_replicates) #specifying singlestages is optional
groups$a -> grouppop_tab
groups$b -> grouppop_plottED
groups$c -> groupweight_plottED
rm(groups)

## Biomass production
production_plottED <- production(output_replicates)

# Save bundle of tabs and plots as RData
save(cleanoutput,
     output_replicates,offspring_replicates,
     vegmass_plottED,repmass_plottED,biomass_tab,
     pop_tab, abund_plottED,
     weekstruct_tab, weekstruct_plottED, relativestruct_tab, relativestruct_plottED, 
     spprichness_tab, spprichness_plottED,
     relabund_tab, rankabund_plottED,
     grouppop_plottED, grouppop_tab, groupweight_plottED,
     file = file.path(analysEDdir, 
                      paste(parentsimID, ".RData", sep = "")))

# Plot all graphs
if(plotall){
  EDplots <- objects(name = environment(), all.names = FALSE, pattern = "_plottED$")
  for(p in seq(1, length(EDplots))){
    ggsave(filename = file.path(analysEDdir, paste(p, ".jpeg", sep = "")) , 
           plot = mget(EDplots)[[p]], #get the actual object (p has the variable name)
           device = "jpeg")}
}