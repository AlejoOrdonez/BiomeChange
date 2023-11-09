rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
require(tidyverse)

#setwd("/Volumes/WDElements20TB/Data/GLOBAL/CLIMATE/FUTURE/CMIP5_Climate_copernicus/Data/Processed/PerModelSummary")
setwd("/Volumes/Crucial X6/Data/GLOBAL/CLIMATE/FUTURE/CMIP5_Climate_copernicus/Data/Processed/PerModelSummary")
OutDir <- "~/Documents/GitHub/BiomeChange/Data/Raw/"

# Create a 50 x 50 km raster
WGSRast <- rast(nrows=180, ncols=360) 
eck4Rast <- project(WGSRast,"+proj=eck4")
eck4Rast <- rast(extent=ext(eck4Rast),resolution=50000,crs="+proj=eck4")
rm(WGSRast);gc()

for(DirUse  in c("Historical", "RCP26", "RCP45", "RCP60", "RCP85")){#(DirUse <- "Historical")
  for(FileUse in dir(paste0("./",DirUse,"/Seasonal"))){#(FileUse <- dir(paste0("./",DirUse,"/Seasonal"))[1])
    RastIn <- rast(paste0("./",DirUse,"/Seasonal/",FileUse))
    ClimMnPrj <- project(RastIn,"+proj=eck4",method="near",threads=50)
    # resample to 50 x 50 km
    ClimMn_50x50 <- lapply(1:dim(ClimMnPrj)[3],
                           function(x){
                             round(resample(ClimMnPrj[[x]],eck4Rast,threads=50),ifelse(x<=4,0,1))
                           })
    ClimMn_50x50 <- do.call("c",ClimMn_50x50)
    writeRaster(ClimMn_50x50,paste0(OutDir,"/By_Model/",DirUse,"/",FileUse),overwrite=TRUE)
    rm(list = c("ClimMn_50x50","ClimMnPrj","RastIn"));gc()
  }
  print(DirUse)
}





