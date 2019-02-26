#library(tidyverse)
#library(rasterVis)
#library(grid)
#library(gtable)

creatEDnlm <- function(loss, size, dir, createcontrol = TRUE){

    # required packages
    require(NLMR)
    require(landscapemetrics)
    require(landscapetools)
    require(raster)

    # create random cluster nlm
    fragmented <- nlm_randomcluster(ncol = size, nrow = size,
                                    p = 0.1,
                                    ai = c(loss, 1-loss))

				    
    # write files
    writeRaster(fragmented,
                paste(file.path(dir,"frag_"),loss*100, "_", size,".grd", sep = ""),
                format = "raster")

if (createcontrol){
control = fragmented
values(control) = 1
# write control file
writeRaster(control,
                paste(file.path(dir,"control_"),loss*100, "_", size,".grd", sep = ""),
                format = "raster")
}

    }
