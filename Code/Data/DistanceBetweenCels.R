# Distances 0.25ArcMin
rm(list=ls());gc()
require(terra)
require(maptools)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
data(wrld_simpl)
wrld_simpl2 <- vect(wrld_simpl)
wrld_simpl2 <- project(wrld_simpl2,"EPSG:4326")
qddRast <- rast(nrows=180*4, ncols=36*4,nlyrs=2,crs="EPSG:4326")

values(qddRast) <- crds(qddRast)
qddRast <- c(qddRast,(qddRast+0.25))

LongDist <- app(qddRast,
         function(x){
           out <- distance(matrix(x[c(1,2)],ncol=2),
                           matrix(x[c(3,2)],ncol=2),
                           lonlat=T)   
           return(out)
           })
plot(LongDist/1000)
quantile(LongDist[]/1000,na.rm=T)
writeRaster(LongDist, "./Data/CMIP5/Processed/CellDist/X.25ArcMIn/LongDist.tif", overwrite = TRUE)

LatDist <- app(qddRast,
         function(x,na.rm=T){
           out <- distance(matrix(x[c(1,2)],ncol=2),
                    matrix(x[c(1,4)],ncol=2),
                    lonlat=T)   
           return(out)
         })
plot(LatDist/1000)
quantile(LatDist[]/1000,na.rm=T)
writeRaster(LatDist, "./Data/CMIP5/Processed/CellDist/X.25ArcMIn/LatDist.tif", overwrite = TRUE)

# Distances 1ArcMin

rm(list=ls());gc()
require(terra)
require(maptools)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
data(wrld_simpl)
wrld_simpl2 <- vect(wrld_simpl)
wrld_simpl2 <- project(wrld_simpl2,"EPSG:4326")
qddRast <- rast(nrows=180, ncols=36,nlyrs=2,crs="EPSG:4326")

values(qddRast) <- crds(qddRast)
qddRast <- c(qddRast,(qddRast+0.25))

LongDist <- app(qddRast,
                function(x){
                  out <- distance(matrix(x[c(1,2)],ncol=2),
                                  matrix(x[c(3,2)],ncol=2),
                                  lonlat=T)   
                  return(out)
                })
plot(LongDist/1000)
quantile(LongDist[]/1000,na.rm=T)
writeRaster(LongDist, "./Data/CMIP5/Processed/CellDist/X.1ArcMIn/LongDist1ArcMin.tif", overwrite = TRUE)

LatDist <- app(qddRast,
               function(x,na.rm=T){
                 out <- distance(matrix(x[c(1,2)],ncol=2),
                                 matrix(x[c(1,4)],ncol=2),
                                 lonlat=T)   
                 return(out)
               })
plot(LatDist/1000)
quantile(LatDist[]/1000,na.rm=T)
writeRaster(LatDist, "./Data/CMIP5/Processed/CellDist/X.1ArcMIn/LatDist1ArcMin.tif", overwrite = TRUE)



