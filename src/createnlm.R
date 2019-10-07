#' Create a random landscape for control and loss/fragmentation experiments.
#'
#' @param area The total area of the simulation arena, in square meters.
#' @param loss A number [0,1], describing the amount of habitat that is unavailable at a given point, i.e., the proportion of cells that should be assigned 0 in the random cluster.
#' @param inputsdir Path to where the files should be written.
#' @return A raster file, \code{fragmented_area}, the size of \code{area}, with \code{loss \times area} 0s, and a raster file, \code{control_area}, the size of \code{area}, with only 1s.
creatEDnlm <- function(area, loss, inputsdir){

    # required packages
    library(NLMR)
    library(landscapemetrics)
    library(landscapetools)
    library(raster)

    # create random cluster nlm
    # -------------------------
    ncells = area # n of cells of 1m² covered by the landscape areaa
    side_length <- as.integer(round(sqrt(ncells)))
    fragmented <- nlm_randomcluster(ncol = side_length, nrow = side_length,
    	       	  	            resolution = 1,
                                    p = 0.1,
                                    ai = c(loss, 1-loss),
				    rescale=TRUE)

				    
    # write files
    # -----------
    writeRaster(fragmented,
                file.path(inputsdir, paste("frag_",loss*100, "_", area,"m2.grd", sep = "")),
                format = "raster")

    # create control file, if there isn't one yet
    controlfilename = paste("control_", area,"m2.grd", sep = "")
    if (!(controlfilename %in% list.files(inputsdir))){
        control = fragmented
        values(control) = 1 #size is the same but values differ
        # write control file
        writeRaster(control,
                    file.path(inputsdir, controlfilename),
                    format = "raster")
    }

}
