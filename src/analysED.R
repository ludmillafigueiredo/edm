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
#library(FactoMineR);
#library(factoextra); #might require mvtnorm 1.0-6, which does not require R 3.5
		     #devtools::install_version("mvtnorm", version = "1.0-6",
		     #                          repos = "http://cran.us.r-project.org")
library(corrplot);
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
allocation_plot <- biomass_allocation(statevars, analysEDMdir)
growthcurve_plot <- growthcurves(statevars, analysEDMdir)

#### Spp. abundances #####
out.list <- pop_abundances(statevars, analysEDMdir)
pop_vars <- out.list$tab
spabund_plot <- out.list$plot

#### Spp. population structure #####
out.list <- pop_structure(pop_vars, seed_production, simID, analysEDMdir)
absltstruct_plot <- out.list$absplot
rltvstruct_plot <- out.list$rltvplot

#### Species richness ####
grouprichness_plot <- spp_richness(statevars, pop_vars, simID, disturbance, analysEDMdir)

#### Species rank-abundance ####
timesteps <- factor(c(min(pop_vars$week), max(pop_vars$week)))
rankabund_plot <- rank_abundances(pop_vars, timesteps, analysEDMdir)

#### Biomass production ####
production_plot <- biomass_production(statevars, analysEDMdir)

#### Life history ####
lifeevents_plot <- lifehistory(outputsdir, simID, analysEDMdir)

#### Age traits ####
out.list <- agetraits(statevars, analysEDMdir)
meanage_tab <- out.list$tab
meanage_jplot <- out.list$jplot
meanage_aplot <- out.list$aplot

#### Seed dynamics ####
#seeddynamics(statevars, seed_production, analysEDMdir)

#### Plot diagnostics ####
diagnostic.grid <- plot_grid(allocation_plot, growthcurve_plot, spabund_plot, absltstruct_plot, rltvstruct_plot, grouprichness_plot, rankabund_plot, production_plot, meanage_jplot, meanage_aplot,
nrow = 5, ncol = 2)
save_plot(file.path(analysEDMdir,"diagnostic.png"), plot = diagnostic.grid,
 nrow = 5, ncol  = 2,
 base_height = 5, base_width = 7,
 dpi = 600)

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
