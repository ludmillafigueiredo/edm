#' Analysis of change in trait values distribution
traitchange  <- function(statevars, timesteps, species){
  
  if(missing(species)){
    species <- unique(statevars$sp)
  }
  species <- factor(species)
  
  #Plot trait distribution for the species, at different time steps
  traitvalue_tab <- statevars%>%
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
  traitsummary_tab  <- statevars%>%
    dplyr::select(-c(id, stage, kernel, clonality, age, xloc, yloc, 
              leaves, stem, root, repr, mated, simID))%>%
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
