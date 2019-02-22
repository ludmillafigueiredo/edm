# !Run it from the model's directory!
# This is an R script intended at faster analysis. Details on the RNotebook file of same name

# Load packages
library(tidyverse);
library(viridis);
library(grid);
library(gridExtra);

# Get outputs
indoutput <- function(input,simID,cluster,outdir){
  
  inid <- file.path(paste(input,".csv", sep = ""))
  outid <- file.path(simID); #simID
  
  # cluster refers to where the analusis is happening, not where the simulation ran
  if (cluster == "gaia") { 
    home = file.path("/home/luf74xx/Dokumente/model")
  } else if (cluster == "hpc") {
    home <- indir <- file.path("/home/ubuntu/model");
  } else {
    home <- file.path("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model");
  }
  
  indir <- file.path(home,"inputs");
  outdir <- file.path(home, outdir);
  
  # Outputs:
  #outdir <- file.path("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model/EDoutputs") #Thinkpad
  outraw <- as_tibble(read.table(file.path(outdir, outid, "orgsweekly.txt"), header = TRUE, sep = "\t"));
  # Clean raw data
  ## Take parentheses out of location column ("()")
  loc <- gsub("[\\(|\\)]", "", outraw$location)
  ## split location in 3 columns
  loc <- as.data.frame(matrix(unlist(str_split(loc, ",")),ncol =3,byrow = T))
  names(loc) = c("xloc", "yloc", "floc")
  ## Complete and clean table
  outdata <- as_tibble(cbind(select(outraw,-one_of("location")),
                             loc))
  write.csv(outdata, file.path(outdir, simID, paste(simID, "idout.csv", sep = "")), row.names = FALSE)
  rm(loc)
  return(list(a = outdata, b = outdir, c = indir))
}
# Extract and plot species-specific biomasses mean and variation for each stage
adultmass <- function(outdata,plotit){
  biomass <- outdata%>%
    select(week,id,stage,sp,veg,repr)
  # Filter all adult individuals and use colours to distinguish the lines representing each individual
  if (plotit == TRUE){
    indvegmass.plot <- ggplot(filter(biomass, stage == "a"), aes(x=week, y= veg, color = factor(id)))+
      geom_line() +
      geom_point()+
      labs(x = "Week", y = "Individual vegetative biomass (g)",
           title = "Adult weekly weight")+
      theme(legend.position = "none")
    
    indrepmass.plot <- ggplot(biomass, aes(x=week, y= repr, color = factor(id)))+
      geom_line() +
      geom_point()+
      labs(x = "Week", y = "Total reproductive biomass (g)",
           title = "Individual reproductive biomass variation")+
      theme(legend.position = "none")
  }
  return(list(a = indvegmass.plot, b = indrepmass.plot, c = biomass))
}
# Extract and plot species-specific biomasses mean and variation for each stage
sppstagemass <- function(biomass, spp){
  
  if (missing(spp)){
    spmass <- biomass%>%
      group_by(week,sp,stage)%>%
      summarise(veg.mean = mean(veg), veg.sd = sd(veg), repr.mean = mean(repr), repr.sd = sd(repr))
  }else{
    spps <- paste("p",spp,sep="-")
    spmass <- biomass%>%
      group_by(week,sp,stage)%>%
      filter(sp %in% spps)%>%
      summarise(veg.mean = mean(veg), veg.sd = sd(veg), repr.mean = mean(repr), repr.sd = sd(repr))
  }
  
  pd <- position_dodge(0.1)
  
  vegmass.plot <- ggplot(spmass, aes(x=week, y= veg.mean, color = factor(stage)))+
    geom_errorbar(aes(ymin=veg.mean-veg.sd, ymax=veg.mean+veg.sd), colour = "black", width=.01, position=pd)+
    geom_line(position=pd)+
    geom_point(position=pd)+
    labs(x = "Week", y = "Mean vegetative biomass (g)",
         title = "Annual biomass variation",
         subtitle = "Mean vegetative biomass (+-sd) per species")+
    facet_wrap(~sp,nrow=5)+
    theme(legend.position = "bottom")
  
  repmass.plot <- ggplot(spmass, aes(x=week, y= repr.mean,  colour = factor(sp)))+
    geom_errorbar(aes(ymin=repr.mean-repr.sd, ymax=repr.mean+repr.sd), 
                  colour = "black", width=.01,position=pd)+
    geom_line()+
    labs(x = "Week", y = "Mean reproductive biomass (g)",
         title = "Annual reproductive biomass variation",
         subtitle = "Mean reproductive biomass (+-sd) per species")+
    theme(legend.position = "bottom")
  
  return(list(a = vegmass.plot, b = repmass.plot, c = spmass))
}
# Extract and plot population abundances
popabund <- function(biomass){
  popdata <- biomass%>%
    select(week,sp,stage) #separated for outputting
  spabund <- popdata%>%
    group_by(week,sp)%>%
    summarize(abundance = n())
  
  abund.plot <- ggplot(spabund, aes(x = week, y = abundance, color = factor(sp)))+
    geom_line() + 
    geom_point() +
    labs(x = "Week", y = "abundanceance (n individuals)",
         title = "Species abundance variation")+
    theme(legend.position = "none")
  
  return(list(a = popdata, b = abund.plot))
}
abund_rep <- function(simID){
  
  orgs <- numeric(0)
  
  for(simID in c(paste(simID, seq(1,reps), sep = "_"))) {
    load(file.path(simID,paste(simID,".RData", sep = "")))
    #include a replication column, to identify replicates
    outdatasim <- mutate(outdata, repli = as.factor(rep(simID, nrow(outdata))))
    if(length(orgs) == 0){
      orgs <- outdatasim
    }else{
      orgs <- bind_rows(orgs, outdatasim, .id = "repli")
    } 
  }
  orgs <- orgs[,-1]
  # Abundance
  popdata <- orgs%>%
    select(week,sp,stage,repli)
  
  spabund <- popdata%>%
    group_by(week,sp,repli)%>%
    summarize(abundance = n())%>%
    ungroup()%>%
    group_by(week,sp)%>%
    summarize(abundance.mean = mean(abundance), abundance.sd = sd(abundance, na.rm = TRUE))
  
  pd <- position_dodge(0.1)
  abundplot <- ggplot(spabund, aes(x = week/52, y = abundance.mean, group = sp, color = factor(sp)))+
    geom_errorbar(aes(ymin = abundance.mean-abundance.sd, ymax= abundance.mean+abundance.sd), 
                  stat = "identity", colour = "black", width=.01, position=pd)+
    geom_line(position=pd)+
    geom_point(position=pd)+
    labs(x = "Year", y = "Abundance (mean +- sd)",
         title = "Species abundance variation")+
    theme(legend.position = "none")
  
  return(list(a = as.data.frame(popdata), b = abundplot))
}
# Extract and plot populations structure
popstruct <- function(popdata,offspringfile,simID,spp){
  
  if (missing(spp)){
    # Population strucuture: abundances
    weekstruct <- popdata%>%
      group_by(week, sp, stage)%>%
      summarize(abundance = n())%>%
      bind_rows(.,select(offspringfile,-mode)) #merge offspring file
    # Population strucuture: proportions
    propstruct.tab <- popdata%>%
      group_by(week,stage)%>%
      summarize(abundance = n())%>%
      complete(week,stage,fill = list(abundance = 0))%>%
      group_by(week)%>%
      mutate(proportion = 100*abundance/sum(abundance))%>%
      complete(week,stage,fill = list(proportion = 0))
  }else{
    spps <- paste("p",spp,sep="-")
    # Population strucuture: abundances
    weekstruct <- popdata%>%
      filter(sp %in% spps)%>%
      group_by(week, sp, stage)%>%
      summarize(abundance = n())
    # Population strucuture: proportions
    propstruct.tab <- popdata%>%
      filter(sp %in% spps)%>%
      group_by(week,stage)%>%
      summarize(abundance = n())%>%
      complete(week,stage,fill = list(abundance = 0))%>%
      group_by(week)%>%
      mutate(proportion = 100*abundance/sum(abundance))%>%
      complete(week,stage,fill = list(proportion = 0))
  }
  # graphs
  weekstruct.plot <- ggplot(weekstruct, aes(x = week, y = abundance, color= factor(stage)))+
    geom_line()+
    geom_point()+
    labs(x = "Abundance", 
         y = "Week",
         title = "Population structure")+
    #scale_color_discrete("Stages:", labels = c("Adults", "Seeds", "Juveniles"))+
    facet_wrap(~sp, nrow = length(unique(weekstruct$sp)))+
    theme(legend.position = "none")+
    theme_minimal()
  
  propstruct.plot <- ggplot(propstruct.tab,
                            aes(x = week, y =proportion, color = as.factor(stage)))+
    geom_line()+
    geom_point()
  
  return(list(a = weekstruct, b = weekstruct.plot,
              c = propstruct.tab, d = propstruct.plot))
}
# Species richness
richness <- function(popdata,simID,disturbance,tdist){
  # extract richness from output
  spprichnesstab <- popdata%>%
    group_by(week)%>%
    summarize(richness = length(unique(sp)))
  
  # plot it. Identify disturbance, when it happens
  if (disturbance == "none" && missing(tdist)){
    spprichnessplot <- ggplot(spprichnesstab, aes(x = week, y = richness))+
      geom_line(color = "dodgerblue2", size = 1.25)+
      geom_point()+
      theme_minimal()
  } else {
    
    if (disturbance == "arealoss10"){
      text <- "10% area lost"
    } else if (disturbance == "arealoss50"){
      text <- "50% area lost"
    } else if (disturbance == "arealoss90") {
      text <- "90% area lost"
    } else {
      text <- "Disturbance"
    }
    
    #my_grob <- grobTree(textGrob(text, x = 0.5, y = 0.9, hjust = 0, 
    #                             gp = gpar(fontsize = 10, fontface = "italic")))
    spprichnessplot <- ggplot(spprichnesstab, aes(x = week, y = richness))+
      geom_line(color = "dodgerblue2", size = 1.25)+
      geom_point()+
      geom_vline(xintercept = tdist, linetype = 2, color = "red")+
      labs(x = "Time", y = "Spp. richness")+
      theme_minimal()
    #annotation_custom(my_grob)+
    
  }
  return(list(a = spprichnesstab, b = spprichnessplot))
}
richness_reps <- function(popdata,simID,disturbance,tdist){
  
  # get reps
  orgs <- numeric(0)
  for(simID in c(paste(simID, seq(1,25), sep = "_"))) {
    load(file.path(simID,paste(simID,".RData", sep = "")))
    #include a replication column, to identify replicates
    outdatasim <- mutate(outdata, repli = as.factor(rep(simID, nrow(outdata))))
    if(length(orgs) == 0){
      orgs <- outdatasim
    }else{
      orgs <- bind_rows(orgs, outdatasim, .id = "repli")
    } 
  }
  orgs <- orgs[,-1]
  popdata <- orgs%>%
    select(week,sp,stage,repli)
  
  # extract richness from output
  spprichnesstab <- popdata%>%
    group_by(week, repli)%>%
    summarize(richness = length(unique(sp)))%>%
    ungroup()%>%
    group_by(week)%>%
    summarize(richness.mean = mean(richness), richness.sd = sd(richness))
  
  # plot it. Identify disturbance, when it happens
  pd <- position_dodge(0.1)
  if (disturbance == "none" && missing(tdist)){
    spprichnessplot <- ggplot(spprichnesstab, aes(x = week/52, y = richness.mean))+
      geom_errorbar(aes(ymin = richness.mean-richness.sd, ymax = richness.mean+richness.sd),
                    color = "black", width=.01, position=pd)+
      geom_line(color = "dodgerblue2", size = 1.25, position = pd)+
      labs(x = "Year", y = "Spp. richness")+
      geom_point(position = pd)+
      #ylim(0,10)+
      theme_minimal()
  } else {
    
    ## Anotate the graph with disturbance type
    # if (disturbance == "arealoss10"){
    #   text <- "10% area lost"
    # } else if (disturbance == "arealoss50"){
    #   text <- "50% area lost"
    # } else if (disturbance == "arealoss90") {
    #   text <- "90% area lost"
    # } else {
    #   text <- "Disturbance"
    # }
    #my_grob <- grobTree(textGrob(text, x = 0.5, y = 0.9, hjust = 0, 
    #                             gp = gpar(fontsize = 10, fontface = "italic")))
    
    spprichnessplot <- ggplot(spprichnesstab, aes(x = week/52, y = richness.mean))+
      geom_line(color = "dodgerblue2", size = 1.25, position = pd)+
      geom_point(position = pd)+
      geom_vline(xintercept = tdist, linetype = 2, color = "red")+
      labs(x = "Time", y = "Spp. richness")+
      ylim(0,10)+
      theme_minimal()
    #annotation_custom(my_grob)+
    
  }
  return(list(a = spprichnesstab, b = spprichnessplot))
}
# Species rank-abundance
rankabund <- function(abund.tab){
  # get relative abundances
  rel_abund.tab <- abund.tab%>%
    ungroup%>%
    filter(week == 448)%>%
    group_by(sp)%>%
    summarize(abundance = n())%>%
    ungroup%>%
    mutate(relabund = abundance/sum(abundance))
  
  rel_abund.plot <- ggplot(rel_abund.tab, aes(x = reorder(sp, -relabund), y = relabund))+
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_text(angle = 50, size = 10, vjust = 0.5))
  
  return(list(a = rel_abund.tab, b = rel_abund.plot))
}
rankabund_reps <- function(abund.tab,t){
  
  # summarized sps abundances 
  
  # get relative abundances
  rel_abund.tab <- abund.tab%>%
    ungroup%>%
    filter(week == t)%>%
    group_by(sp, repli)%>%
    summarize(abundancer = n())%>%
    ungroup()%>%
    group_by(sp)%>%
    summarize(abundance = mean(abundancer))%>%
    ungroup%>%
    mutate(relabund = abundance/sum(abundance))
  
  rel_abund.plot <- ggplot(rel_abund.tab, aes(x = reorder(sp, -relabund), y = relabund))+
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_text(angle = 50, size = 10, vjust = 0.5))
  
  return(list(a = rel_abund.tab, b = rel_abund.plot))
}
# Biomass production
production <- function(biomass, structured, area){
  prod <- biomass%>%
    select(week,veg,repr)%>%
    group_by(week)%>%
    summarize(prod.ha = sum(veg,repr)/1000000/area) #g to T and /1ha conversions
  
  total <- ggplot(prod, aes(x= week, y = prod.ha))+
    geom_point()+
    geom_line()+
    labs(x = "Week",
         y = "Tones/ha")+
    theme_minimal()
  
  if (structured != TRUE){
    require(gridExtra)
    structprod <- biomass%>%
      select(week,stage,veg,repr)%>%
      group_by(week,stage)%>%
      summarize(prod.ha = sum(veg,repr)/1000000/area)
    
    stagesprod <- ggplot(structprod, aes(x= week, y = prod.ha, color = factor(stage)))+
      geom_point()+
      geom_line()+
      labs(x = "Week",
           y = "Tones/ha")+
      theme_minimal()
    
    return(list(a = total, b = stagesprod))
    
  }else{
    return(total)
  }
}
# Population structure by group size
group.pop <- function(outdata, offspringfile, singlestages){
  
  # create table
  grouppop.tab <- outdata %>%
    group_by(week, sp,e_mu,stage)%>%
    summarize(abundance = n())%>%
    ungroup()%>%
    bind_rows(.,select(offspringfile,-mode))%>% # merge seed info
    group_by(sp)%>% # necessary to fill in seed size according to species
    fill(e_mu) 
  
  #create plot
  if (missing(singlestages)){
    grouppop.plot <- ggplot(data = grouppop.tab,
                            aes(x = week, y = abundance, colour = stage))+
      geom_line()+
      geom_point()+
      facet_wrap(~e_mu, ncol = 1, nrow = 3)+
      labs(title = "Population structure per group size")+
      theme_bw()
  }else{
    grouppop.plot <- ggplot(data = grouppop.tab %>%
                              filter(stage %in% singlestages),
                            aes(x = week, y = abundance, colour = stage))+
      geom_line()+
      geom_point()+
      facet_wrap(~e_mu, ncol = 1, nrow = 3)+
      labs(title = "Population structure per group size")+
      theme_bw()
  }
  return(list(a = grouppop.plot, b = grouppop.tab))
}

# Set up directory to store analysis
dir.create(file.path(getwd(),paste(simID,"all", reps,sep = "_")))

outputs <- indoutput(input,simID,cluster,outdir)
outputs$a -> outdata  
outputs$b -> indir
rm(outputs)

# Individual vegetative and reproductive biomasses
mass <- adultmass(outdata,TRUE)
mass$a -> indvegmass.plot 
mass$b -> indrepmass.plot 
mass$c -> biomass.tab
rm(mass)

# Stage biomass variation, summarized per stage, for each species
stgmass <- sppstagemass(biomass.tab)
stgmass$a -> vegmass.plot
stgmass$b -> repmass.plot
stgmass$c -> stgmasssummary.tab
rm(stgmass)

# Species abundance variation
abund <- popabund(biomass.tab)
abund$a -> abund.tab
abund$b -> abund.plot
rm(abund)

# Population structure variation
pop <- popstruct(abund.tab,offspringfile,simID)
pop$a -> popstruct.tab
pop$b -> popstruct.plot
rm(pop)

# Species richness
rich <- richness(abund.tab,simID,disturbance)
rich$a -> spprichness.tab 
rich$b -> spprichness.plot
rm(rich)

# Species rank-abundance 
rank <- rankabund(abund.tab)
rank$a -> rank.tab
rank$b -> rank.plot
rm(rank)

# Biomass production/ha
prod <- production(biomass.tab, structured = TRUE, area = 5)

# Population structure by group size
group <- group.pop(outdata, offspringfile)
group$a -> grouppop.plot
group$b -> grouppop.tab
rm(group)

# Weight by group size 
group.weight <- ggplot(data = outdata %>%
                         group_by(week,e_mu,stage)%>%
                         summarize(weight = mean(veg), weight.sd = sd(veg)),
                       aes(x = week, y = weight, colour = stage))+
  geom_errorbar(aes(ymin = weight-weight.sd, ymax = weight + weight.sd), 
                colour = "black", width=.01, position = position_dodge(0.1))+
  geom_line(position = position_dodge(0.1))+
  geom_point(position = position_dodge(0.1))+
  facet_wrap(~e_mu, ncol = 1, nrow = 3)+
  theme_bw()

## SAVE BUNDLE OF GRAPHS AND TABLES AS RDATA
save(outdata,
     indvegmass.plot,indrepmass.plot,
     vegmass.plot,repmass.plot,biomass.tab,stgmasssummary.tab,
     abund.plot,abund.tab,
     popstruct.tab,popstruct.plot,
     spprichness.tab, spprichness.plot,
     prod,
     grouppop.plot,
     grouppop.tab,
     group.weight,
     rank.tab,
     rank.plot,
     file = file.path(simID,paste(simID,".RData", sep = "")))