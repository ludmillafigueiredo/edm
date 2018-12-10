# !Run it from the model's directory!
# This is an R script intended at faster analysis. Details on the RNotebook file of same name

# Load packages
library(tidyverse);
library(viridis);
library(grid);
library(gridExtra);

# Get outputs
indoutput <- function(input,simID,cluster){
  
  inid <- file.path(paste(input,".csv", sep = ""))
  outid <- file.path(simID); #simID
  
  if (cluster == "gaia") {
    home = file.path("/home/luf74xx/Dokumente/model")
  } else if (cluster == "hpc") {
    home <- indir <- file.path("/home/ubuntu/model");
  } else {
    home <- file.path("/home/ludmilla/Documents/uni_wuerzburg/phd_project/thesis/model");
  }
  
  indir <- file.path(home,"inputs");
  outdir <- file.path(home,"EDoutputs");
  
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
  write.csv(outdata, file.path(outdir, simID,paste(simID, "idout.csv", sep = "")), row.names = FALSE)
  rm(loc)
    return(list(a = outdata, b = outdir, c = indir))
}
# Extract and plot individual vegetative and reproductive biomasses.
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
    labs(x = "Week", y = "Abundance (n individuals)",
         title = "Species abundance variation")+
    theme(legend.position = "none")
  
  return(list(a = popdata, b = abund.plot))
}
# Extract and plot populations structure
popstruct <- function(popdata,simID,spp){
  
  if (missing(spp)){
    # Population strucuture: abundances
    weekstruct <- popdata%>%
      group_by(week, sp, stage)%>%
      summarize(abundance = n())
    # Population strucuture: proportions
    propstruct.tab <- popdata%>%
      group_by(week,stage)%>%
      summarize(abund = n())%>%
      complete(week,stage,fill = list(abund = 0))%>%
      group_by(week)%>%
      mutate(proportion = 100*abund/sum(abund))%>%
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
      summarize(abund = n())%>%
      complete(week,stage,fill = list(abund = 0))%>%
      group_by(week)%>%
      mutate(proportion = 100*abund/sum(abund))%>%
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
    distinct(week,sp)%>%
    summarize(richness = n())
  
  # plot it. Identify disturbance, when it happens
  if (disturbance == "none" && missing(tdist)){
    spprichnessplot <- ggplot(spprichnesstab, aes(x = week, y = richness))+
      geom_line(color = "dodgerblue2", size = 1.25)+
      geom_point()+
      geom_vline(xintercept = 50, linetype = 2, color = "red")+
      theme_minimal()
  } else {
    
    if (disturbance == "arealoss10"){
      text <- "10% area lost"
    } else if (disturbance == "arealoss50"){
      text <- "50% area lost"
    } else if (disturbance == "arealoss90") {
      text <- "90% area lost"
    } 
    
    #my_grob <- grobTree(textGrob(text, x = 0.5, y = 0.9, hjust = 0, 
    #                             gp = gpar(fontsize = 10, fontface = "italic")))
    
    spprichnessplot <- ggplot(spprichnesstab, aes(x = week, y = richness))+
      geom_line(color = "dodgerblue2", size = 1.25)+
      geom_point+
      geom_vline(xintercept = tdist, linetype = 2, color = "red")+
      labs(x = "Time", y = "Spp. richness")+
      theme_minimal()
      #annotation_custom(my_grob)+
     
  }
  return(list(a = spprichnesstab, b = spprichnessplot))
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

# Set up directory to store analysis
dir.create(file.path(getwd(),simID))

outputs <- indoutput(input,simID,cluster)
outputs$a -> outdata  
outputs$b -> outdir
outputs$c -> indir
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

# Species structure variation
pop <- popstruct(abund.tab,simID)
pop$a -> popstruct.tab
pop$b -> popstruct.plot

# Species richness
rich <- richness(abund.tab,simID,disturbance)
rich$a -> spprichness.tab 
rich$b -> spprichness.plot
rm(rich)

# Biomass production/ha
prod <- production(biomass.tab, structured = TRUE, area = 5)

## SAVE BUNDLE OF GRAPHS AND TABLES AS RDATA
save(outdata,
     indvegmass.plot,indrepmass.plot,
     vegmass.plot,repmass.plot,biomass.tab,stgmasssummary.tab,
     abund.plot,abund.tab,
     popstruct.tab,popstruct.plot,
     spprichness.tab, spprichness.plot,
     prod,
     file = file.path(simID,paste(simID,".RData", sep = "")))

# Save graphic outputs 
# ggsave("richness.eps", spprichness.plot, 
#        device = "eps", getwd(),
#        width = 15, height = 10, units = "cm")













