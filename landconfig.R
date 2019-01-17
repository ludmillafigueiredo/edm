require(dplyr)
require(tidyr)
require(gstat)
require(rgdal)
require(raster)
require(rgeos)
require(scales)
require(stringr)
options(stringsAsFactors=FALSE)

# read in shape file of initial landscape
initialshape <- readOGR(file.path(initialland))
initialdistances <- pointDistance(initialshape, lonlat = FALSE) # lonlat is TRUE if coordinates are in degrees
initialareas <- area(initialshape)

# read in shape file of disturbed landscape
disturbshape <- readOGR(file.path(disturbland))
disturbdistances <- pointDistance(disturbshape, lonlat = FALSE) # lonlat is TRUE if coordinates are in degrees
disturbareas <- area(disturbshape)

# find middle point of centroids and create a buffer around it
