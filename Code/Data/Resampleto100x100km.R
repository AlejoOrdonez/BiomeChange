require(terra)
e <- rast("/Volumes/Crucial X6/Data/GLOBAL/CLIMATE/FUTURE/CMIP5_Climate_copernicus/Data/Processed/AllModelSumm1Arcmin/RCP26/BIOCLIM/AllModels_RCP26_1980.tif")
 plot(e[[1]])
a <- rast(nrows=180, ncols=360)
b<-project(e[[1]],"+proj=eck4")
plot(b)
d <- rast(extent=ext(b),resolution=100000)

f <- resample(b,d)
plot(f)
plot(e[[1]])


e <- rast(nrows=180, ncols=360) 
b<-project(e[[1]],"+proj=eck4")
d <- rast(extent=ext(b),resolution=100000,crs="+proj=eck4")
res(59455.05)
