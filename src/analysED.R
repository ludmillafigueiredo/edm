## This script provides summary analysis of the output of a simulation (one replicate):
## It requires
## `simID`: string of complete simulation id (included replicate number)
## `disturbance`: string "none", "area_loss", "poll_loss", etc (see 'inputs.jl' for all options. 
## `tdist`: vector of disturbance time or NULL (`nothing` in Julia)
## `outputsdir`: path to directory where outputs are store
## These are provided by simulation.

#### Load packages ####
library(tidyverse);
library(viridis);
library(readr);
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
analysEDMdir <- file.path(outputsdir, "analysED")
dir.create(analysEDMdir)

#### Curstomized theme, labels, etc ####
source(file.path(EDMdir,"theme_edm.R"))
theme_set(theme_edm())
size_labels = c('0.0001' = "Small",
                '0.0003' = "Medium",
                '0.001' = "Big")
options(scipen = 999)

#### Summary analysis ####
#-------------------------#
source(file.path(EDMdir, "summaryanalysis_functions.R"))
source(file.path(EDMdir, "traitspace_functions.R"))

#### Clean main output ####
statevars <- format_statevars(outputsdir, simID, analysEDMdir)
seed_production <- format_seedprod(outputsdir, simID, analysEDMdir)

#### Growth #####
biomass_allocation(statevars, analysEDMdir)
growthcurves(statevars, analysEDMdir)

#### Spp. abundances #####
pop_vars <- pop_abundances(statevars, analysEDMdir)

#### Spp. population structure #####
pop_structure(pop_vars, seed_production, simID, analysEDMdir)

#### Species richness ####
spp_richness(statevars, pop_vars, simID, disturbance, analysEDMdir)

#### Species rank-abundance ####
timesteps <- factor(c(min(pop_vars$week), max(pop_vars$week)))
rank_abundances(pop_vars, timesteps, analysEDMdir)

#### Biomass production ####
biomass_production(statevars, analysEDMdir)

#### Life history ####
lifehistory(outputsdir, simID, analysEDMdir)

#### Age traits ####
agetraits(statevars, analysEDMdir)

#### Seed dynamics ####
#seeddynamics(statevars, seed_production, analysEDMdir)

#### Trait change ####
### trait values
#traitschange <- traitchange(statevars, timesteps)
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

#### Trait analysis ####

#traitvaluesdir <- file.path(analysEDMdir, "traitvalues")
#dir.create(traitvaluesdir)

#for(sp in 1:length(unique(statevars$sp))){
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
#                            paste(unique(statevars$sp)[sp],"_dist.png", sep = "")),
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
#                            paste(unique(statevars$sp)[sp],"_ts.png", sep = "")),
#              plot = ts,
#              dpi = 300,
#              height = 30,
#              width = 40,
#              units = "cm")
#  }
#}
