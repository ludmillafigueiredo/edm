library(NLMR)
library(landscapemetrics)
library(landscapetools)
library(raster)
#library(tidyverse)
#library(rasterVis)
#library(grid)
#library(gtable)

# 25% area loss
frag25 <- nlm_randomcluster(ncol = 3536, nrow = 3536,
                                     p = 0.1,
                            ai = c(0.25, 0.75))
# 75% area loss
frag50 <- nlm_randomcluster(ncol = 3536, nrow = 3536,
                                     p = 0.1,
                            ai = c(0.5, 0.5))
# 75% area loss
frag75 <- nlm_randomcluster(ncol = 3536, nrow = 3536,
                                     p = 0.1,
                            ai = c(0.75, 0.25))

                                        # write files
writeRaster(frag25,
            "/home/ubuntu/model/inputs/frag25.grd",
            format = "raster")

writeRaster(frag50,
            "/home/ubuntu/model/inputs/frag50.grd",
            format = "raster")

writeRaster(frag75,
            "/home/ubuntu/model/inputs/frag75.grd",
            format = "raster")
