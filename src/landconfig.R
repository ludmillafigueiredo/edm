require(sp)
require(dplyr)
require(tidyr)
require(gstat)
require(rgdal)
require(raster)
require(rgeos)
require(scales)
require(stringr)
options(stringsAsFactors=FALSE)

# 1. Read in shape file of initial landscape
initialshape <- readOGR(file.path(initialland))
# 2. Read in shape file of disturbed landscape
if (!is.null(disturbland)){
    disturbshape <- readOGR(file.path(disturbland))
}else{
    disturbshape = NULL
}
# 3. Read in buffer file
if (!is.null(landbuffer)){
    landbuffer <- readOGR(file.path(landbuffer))
}else{
    landbuffer <- NULL
}
# 4.Store relevant measures
landconfig <-list (
    length(initialshape), # initial number of patches
    area(initialshape), # their areas
    sum(area(initialshape)),
    pointDistance(initialshape, lonlat = FALSE), # distances between their middle points
    ifelse(is.null(disturbshape), NULL, length(disturbshape)), # number of fragments
    ifelse(is.null(disturbshape), NULL, area(disturbshape)), # their areas
    ifelse(is.null(disturbshape), NULL, sum(area(disturbshape))),
    ifelse(is.null(disturbshape), NULL, pointDistance(disturbshape, lonlat = FALSE)), # distances between them
    ifelse(is.null(landbuffer), NULL, area(landbuffer))) # area of buffer including all patches/fragments)


# find middle point of centroids and measure the area of the  buffer around it - descarted, because it can vary, depending on the case. Using a buffer shape, instead.
#mean(coordinates(initialshape)[,1]) # mean of centroids x
#mean(coordinates(initialshape)[,2])# mean of centroids y

rm(initialshape,disturbshape,landbuffer)

