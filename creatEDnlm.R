#library(tidyverse)
#library(rasterVis)
#library(grid)
#library(gtable)

creatEDnlm <- function(loss, size, dir){

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
                paste(file.path(dir,"frag_"),loss, "_", size,".grd", sep = ""),
                format = "raster")
    }
