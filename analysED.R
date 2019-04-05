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
library(FactoMineR);
library(factorextra);
library(corrplot);
library(gganimate); #sudo apt-get install cargo install.packages("gifski")
library(cowplot);

# Default values for analysis
    
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
orgreplicates <- function(parentsimID, nreps){
    
    output_repli <- data.frame()
    offspring_repli <- data.frame()
    
    for(sim in c(paste(parentsimID, seq(1, nreps), sep = "_"))) {
        repfolder <- file.path(outdir, parentsimID, sim); # replicate folder
        # read outputs of juv/adults and of offspring
        outdatasim <- read.csv(file.path(repfolder, paste(parentsimID, "indout.csv", sep = "")), header = TRUE)%>%
            mutate(repli = as.factor(rep(sim, nrow(.)))) #include a replication column, to identify replicates
        offspringsim <- read.table(file.path(repfolder, "offspringproduction.csv"), header = TRUE, sep = "\t")%>%
            mutate(repli = as.factor(rep(sim, nrow(.))))
        
        # fill in cleanoutput object (necessary step because of bind_rows)
        if(length(output_repli) == 0){
            output_repli <- outdatasim
            offspring_repli <- offspringsim
        }else{
                                        # identify the lines corresponding to it with its replicate ID
            output_repli <- bind_rows(output_repli, outdatasim, .id = "repli")
            offspring_repli <- bind_rows(offspringsim, offspringsim, .id = "repli")
            offspring_repli$week = factor(offspring_repli$week)
        } 
    }
    return(list(a = output_repli, b = offspring_repli))
}

#' Extract and plot species-specific mean and sd biomass compartiments (these are means of each simulation' means   
#' @param outdatalist
#' @param plotit Boolean specifying whether graph shouuld be plotted or not
stagemass <- function(output_repli, stages = factor(c("a", "j", "s"))){
    
                                        # get biomasses
    biomass_tab <- output_repli%>%
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
popabund <- function(output_repli){
    
                                        # extract relevant variables
    pop_tab <- output_repli%>%
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
popstruct <- function(pop_tab, offspring_repli, parentsimID, spp){
    
    if (missing(spp)){
                                        # Population strucuture: abundances
        weekstruct_tab <- pop_tab%>%
            group_by(week, sp, stage, repli)%>%
            summarize(abundance = n())%>%
            ungroup()%>%
            bind_rows(., select(offspring_repli, -mode))%>% #merge offspring file
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
            ungroup()%>%
            bind_rows(., select(offspring_repli, -mode))%>% #merge offspring file
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
            geom_errorbar(aes(ymin = mean_richness - sd_richness,
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
rankabund <- function(pop_tab, timesteps){

    relabund_tab <- pop_tab%>%
                filter(week %in% timesteps)%>%
                group_by(week, sp, repli)%>%
                summarize(abundance = n())%>%
                ungroup()%>%
                group_by(sp, week)%>%
                summarize(mean_abundance = mean(abundance))%>%
                ungroup()%>%
		group_by(week)%>%
                mutate(relabund = mean_abundance/sum(mean_abundance))

    plot_rankabund <- function(timestep){
    	ggplot(filter(relabund_tab, week %in% timestep),
		aes(x = reorder(sp, -relabund), y = relabund))+
                geom_bar(stat = "identity")+
                theme_minimal()+
                theme(axis.text.x = element_text(angle = 50, size = 10, vjust = 0.5))
    }

    rankabunds_plot <- timesteps%>%
    	map(. %>% plot_rankabund)
      
    return(list(a = relabund_tab, b = rankabunds_plot))
}

#' Population structure by group size
groupdyn <- function(output_repli, singlestages){
    
    # create table
    grouppop_tab <- output_repli%>%
        group_by(week, sp, emass, stage, repli)%>%
        summarize(abundance = n())%>%
        ungroup()%>%
        bind_rows(., select(offspring_repli, -mode))%>% # merge seed info
        group_by(sp)%>% # necessary to fill in seed size according to species
        fill(emass)%>%
        ungroup()%>%
        group_by(week, emass, stage)%>% 
        summarize(mean_abundance = mean(abundance),
                  sd_abundance = sd(abundance))%>%
        ungroup()
    
    # create plots 
    ## filter data, if necessary
    if (missing(singlestages)){
        popdata_to_plot <- grouppop_tab
        weightdata_to_plot <- output_repli
    }else{
        popdata_to_plot <- grouppop_tab %>%
            filter(stage %in% singlestages)
        weightdata_to_plot <- output_repli%>%
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
        facet_wrap(~emass, ncol = 1, nrow = 3)+
        labs(title = "Population structure per group size")+
        theme_minimal()
    
                                        # plot weigh variation
    groupweight_plot <- ggplot(data = weightdata_to_plot%>%
                                      group_by(week, emass, stage, repli)%>%
                                      summarize(repli_mean_weight = mean(veg))%>%
                                      ungroup()%>%
                                      group_by(week, emass, stage)%>%
                                      summarize(mean_weight = mean(repli_mean_weight), 
                                                sd_weight = sd(repli_mean_weight)),
                                  aes(x = week, y = mean_weight, colour = stage))+
        geom_errorbar(aes(ymin = mean_weight-sd_weight, ymax = mean_weight + sd_weight), 
                      colour = "black", width=.01, position = position_dodge(0.1))+
        geom_line(position = position_dodge(0.1))+
        geom_point(position = position_dodge(0.1))+
        facet_wrap(~emass, ncol = 1, nrow = 3)+
        theme_bw()
    
    return(list(a = grouppop_tab, b = grouppop_plot, c = groupweight_plot))
    
}

#' Calculate biomass production
production <- function(output_repli){
    
    production_tab <- output_repli%>%
        select(week,veg,repr,repli)%>%
        group_by(week, repli)%>%
        summarize(production = sum(veg,repr)/(10^3))%>%
        ungroup()%>%
        group_by(week)%>%
        summarize(mean_prod = mean(production),
                  sd_prod = sd(production))%>%
        ungroup()
    
    production_plot<- ggplot(production_tab, aes(x= week, y = mean_prod))+
        geom_errorbar(aes(ymin = mean_prod - sd_prod,
                          ymax = mean_prod + sd_prod),
                      stat = "identity", position = position_dodge(0.01), colour = "gray50", width = 0.01)+
        geom_point()+
        geom_line()+
        labs(x = "Week",
             y = "Biomass production (kg)")+
        theme_minimal()
    
    return(production_plot)
}

#' Analysis of change in trait values distribution

traitchange  <- function(output_repli, timesteps, species){

    # take the mean values for each individual (they are replicated in each experiment) and then plot the violin plots
    traitvalues_tab  <- output_repli%>%
        select(-c(kernel, clonal, age, xloc, yloc, veg, repr, mated, repli))%>%
        group_by(week, id, stage, sp)%>%
        summarize_all(funs(mean = mean,
                           sd = sd))%>%
        ungroup()

    if(missing(species)){
        species <- unique(traitvalues_tab$sp)
    }
    species <- factor(species)
    
    # internal function for different species and the same time step (map() goes through each spp, for the same timesteps). Internal function because maps needs single function and because column selection inside gather() does not work otherwise, probably piped objects get mixed (mapped spp vector and piped trait table).

    traittab4plot <- gather(traitvalues_tab%>%
                      select(-tidyselect::ends_with("_sd"))%>%
                      filter(week %in% timesteps),
                      key = trait,
                      value = value,
                      emass_mean:maxmass_mean,
                      factor_key = TRUE)
		      
    plottrait <- function(spp){
        ggplot(traittab4plot%>%filter(sp %in% spp),
               aes(x = sp, y = value))+
            geom_violin()+
            #geom_dotplot(color = "grey31", fill = "white", alpha = 0.8)+
            geom_boxplot(width = 0.2)+
            facet_wrap(c("week", "trait"), #facet_grid cannot free y axis
                       nrow = length(unique(timesteps)),
                       ncol = length(select(traitvalues_tab, tidyselect::ends_with("_mean"))),
                       scales = "free_y")+
            background_grid(major = "xy", minor = "none")
}

    traitvalues_plots <- species%>%
        map(. %>% plottrait) 

    #gganimation
    #animate(traitchange_plot +
    #        transition_states(week, transition_length = 2, state_length= 3)+
    #        ease_aes('cubic-in-out'),
    #        width = 800, height = 2000, nframes = 100,
    #        device = "png",
    #        resolution = 300,
    #        renderer = file_renderer(analysEDdir, prefix = paste("animtrait", species, sep = "_"), overwrite = TRUE))
        
    return(list(a = traitvalues_tab, b = traitvalues_plots))
}

#' Analysis of change in trait space

traitspacechange  <- function(traitvalues_tab, timesteps){

    # PCA (internal function because it needs to passed to map)
    traitPCA <- function(timestep){

        traitpca <- PCA(traitvalues_tab%>%
			select(-sp, -stage, -ends_with("_sd"))%>%
            		filter(week %in% timestep)%>%
			select(-week, -id),
                    scale.unit = TRUE,
                    ncp = 5,
                    graph = FALSE)

        return(traitpca)
    }

    #separate PCAs for t0 and tend
    traitpcas <- timesteps%>%
        map(., traitPCA)

    # PCA comparing t0 and end
    timepca <- traitPCA(timesteps)
    timepca_plot <- fviz_pca_biplot(timepca,
                          geom.ind = "point", # show points only (but not "text")
                          pointshape = 21,
                          pointsize = 2.5,
                          fill.ind = factor(select(filter(traitvalues_tab, week %in% timesteps), week)), # color by time
                          col.ind = "black",
                                        #addEllipses = TRUE, # Concentration ellipse,
                          col.var = factor(c("size", "reprd", "reprd", "reprd", "reprd", "reprd", "span", "metab", "metab", "metab","reprd", "size")),
                          repel = TRUE,
                                       legend.title = list(fill = "Time-step", color = "Traits"))
    
    return(list(a = traitpcas, b = timepca, c = timepca_plot))

}


############################################################################
#                          Organize analysis output                        #
############################################################################

cleanoutput <- getoutput(parentsimID, nreps, outdir = outdir, EDdir = EDdir)  

## Identify replicates
replicates <- orgreplicates(parentsimID,nreps)
replicates$a -> output_repli
replicates$b -> offspring_repli
rm(replicates)

## Individual vegetative and reproductive biomasses of juveniles and adults (NOT seeds)
mass <- stagemass(output_repli,TRUE) # spp is an optional argument and stage is set to default
mass$a -> vegmass_plot 
mass$b -> repmass_plot 
mass$c -> biomass_tab
rm(mass)

## Species abundance variation
abund <- popabund(output_repli)
abund$a -> pop_tab
abund$b -> spabund_tab
abund$c -> abund_plot
rm(abund)

## Population structure variation
pop <- popstruct(pop_tab, offspring_repli, parentsimID) #specifying species is optional
pop$a -> weekstruct_tab
pop$b -> weekstruct_plot
pop$c -> relativestruct_tab
pop$d -> relativestruct_plot
rm(pop)

## Species richness
rich <- richness(pop_tab, parentsimID, disturbance) # tdist is optional
rich$a -> spprichness_tab 
rich$b -> spprichness_plot
rm(rich)

## Set up time-steps for which to output derived analysis
#if(!("timesteps" %in% ls())){timesteps <- c(min(pop_tab$week), max(pop_tab$week))}
timesteps <- factor(c(min(pop_tab$week), max(pop_tab$week))) #provided or not, must be converted

## Species rank-abundance 
rank <- rankabund(pop_tab, timesteps)
rank$a -> relabund_tab
rank$b -> rankabunds_plot

rm(rank)

## Population structure by group size
groups <- groupdyn(output_repli) #specifying singlestages is optional
groups$a -> grouppop_tab
groups$b -> grouppop_plot
groups$c -> groupweight_plot
rm(groups)

## Biomass production
production_plot <- production(output_repli)

## Trait change
### trait values
traitschange <- traitchange(output_repli, timesteps)
traitschange$a -> traitvalues_tab
traitschange$b -> traitvalues_plot # no 's' so it can be detected by `plotall`
rm(traitschange)
### trait space
traitspace <- traitspacechange(traitvalues_tab, timesteps)
traitspace$a -> traitpcas
traitspace$b -> timepca
traitspace$c -> timepca_plot

# Save bundle of tabs and plots as RData
save(cleanoutput,
     output_repli,offspring_repli,
     vegmass_plot,repmass_plot,biomass_tab,
     pop_tab, abund_plot,
     # weekstruct_tab, weekstruct_plot, relativestruct_tab, relativestruct_plot, 
     spprichness_tab, spprichness_plot,
     relabund_tab, rankabunds_plot,
     grouppop_plot, grouppop_tab, groupweight_plot,
     traitvalues_tab, traitvalues_plot,
     traitpcas, timepca, timepca_plot
     file = file.path(analysEDdir, 
                      paste(parentsimID, ".RData", sep = "")))

# Plot all graphs
EDplots <- objects(name = environment(), all.names = FALSE, pattern = "_plot$")
map(EDplots,
    ~ ggsave(file.path(analysEDdir, paste(.x, ".png", sep ="")), get(.x)))
