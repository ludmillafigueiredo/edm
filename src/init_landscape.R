require(sp)
require(raster)

# Read raster of initial landscape as matrix (Julia processes matrices, not raster)
initial_matrix <- raster::as.matrix(raster(file.path(initial_land)))
