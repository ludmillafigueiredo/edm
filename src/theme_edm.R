options(ggplot2.continuous.colour = "viridis")
scale_colour_discrete <- scale_colour_viridis_d
update_geom_defaults("point", list(size = 0.3))
update_geom_defaults("errorbar", list(colour = "azure"))
update_geom_defaults("errorbar", list(width = 0.01))

theme_edm <- function(){
    theme_minimal(base_size = 8)%+replace%
            theme(
              axis.text = element_text(size = rel(1)),
	      axis.title.x = element_text(size = rel(1.25)),
	      axis.title.y = element_text(size = rel(1.25),
	      		     		  angle = 90),
	      legend.position="bottom",
	      legend.background = element_rect(fill="transparent", colour=NA),
	      legend.key = element_rect(fill="transparent", colour=NA)
	     )
}