#library(tidyverse)
#library(rasterVis)
#library(grid)
#library(gtable)

#' Create nlm landscape for control and loss/fragmentation experiment. Area in ha
creatEDnlm <- function(loss, area, inputsdir){

    # required packages
    require(NLMR)
    require(landscapemetrics)
    require(landscapetools)
    require(raster)

    # create random cluster nlm
    ncells = area*(10^8)/10000 # n of cells of 10000 cm² covered by the landscape area (ha)
    fragmented <- nlm_randomcluster(ncol = round(sqrt(ncells)), nrow = round(sqrt(ncells)),
                                    p = 0.1,
                                    ai = c(loss, 1-loss))

				    
    # write files: name has format frag/control_percentageloss_initialarea.grd
    writeRaster(fragmented,
                file.path(inputsdir, paste("frag_",loss*100, "_", area,"ha.grd", sep = "")),
                format = "raster")

    # create control file, if there isn't one yet
    controlfilename = paste("control_", area,"ha.grd", sep = "")
    if (!(controlfilename %in% list.files(inputsdir))){
        control = fragmented
        values(control) = 1 #size is the same but values differ
        # write control file
        writeRaster(control,
                    file.path(inputsdir, controlfilename),
                    format = "raster")
    }

}
