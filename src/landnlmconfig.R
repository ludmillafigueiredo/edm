require(sp)
require(raster)
require(NLMR)
require(landscapetools)
require(landscapemetrics)

# 1. Read raster of initial landscape as matrix (Julia processes matrices, not raster)
initialmatrix <- raster::as.matrix(raster(file.path(initialland)))
cat("initialand is ", typeof(initialmatrix), "\n")

# 2. Same for disturbed landscape, if that's the case
disturbmatrix <- "notfrag"
if(typeof(disturbland) == "list" & "disturbland" %in% names(disturbland)){
    disturbmatrix <- raster::as.matrix(raster(file.path(disturbland$disturbland)))
}
#cat("disturbland$disturbland is ", disturbland$disturbland, "\n")
cat("disturbland is ", typeof(disturbland), "\n")
cat("disturbmatrix is ", disturbmatrix, "\n")
