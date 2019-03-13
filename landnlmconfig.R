library(raster)
library(NLMR)
library(landscapetools)
library(landscapemetrics)

# 1. Read raster of initial landscape as matrix (Julia processes matrices, not raster)
initialmatrix <- raster::as.matrix(raster(file.path(initialland)))
# 2. Same for disturbed landscape, if that's the case
if(typeof(disturbland) == "character" ){
    disturbmatrix <- raster::as.matrix(raster(file.path(disturbland)))
}else{
    #if the column does not contain a file path, it either simulates habitat loss, or no land disturbance at all. In both cases, assignement is handled at read_landpars().
disturbmatrix = "none"
}

rm(initialland, disturbland)
# 4. Store everything in the list Julia will rcopy to the model environment
landmatrices <- list(
    initialmatrix, # initial landscape matrix
    disturbmatrix) # disturbed landscape matrix
