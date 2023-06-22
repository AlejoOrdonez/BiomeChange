rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")

# BioclimVars to use
BioclimVars <- c(5, # Max Temperature of Warmest Month
                 6, # Min Temperature of Coldest Month
                 13,# Precipitation of Wettest Month
                 14) # Precipitation of Driest Month

# Create a 100 x100 km raster
WGSRast <- rast(nrows=180, ncols=360) 
eck4Rast <- project(WGSRast,"+proj=eck4")
eck4Rast <- rast(extent=ext(eck4Rast),resolution=50000,crs="+proj=eck4")
rm(WGSRast);gc()

# Load the biome raster
BiomeBsLn <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_50km.tif")

# Estimate the baseline mean and SD of climate variables
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[2])
  # Load the 1980 to 2010 Historical data for an RCP to Create a climate normal raster for each evaluated variable  
  ClimNorm4RCPTmp <- lapply(1980:2010,
                            function(YearUse){#(YearUse <- c(1980:2010)[1])
                              tmp <- rast(paste0("./Data/CMIP5/Processed/AllModelSumm1Arcmin/",RCP,"/BIOCLIM/AllModels_",RCP,"_",YearUse,".tif"),
                                          lyrs= BioclimVars)
                              tmp <- project(tmp,"+proj=eck4")
                              resample(tmp,
                                       eck4Rast,
                                       method = "near")
                            })
  # Estimate the mean of each band for the Climate Normal
  ClimNormMn <- do.call("c",
                        lapply(names(ClimNorm4RCPTmp[[1]]),
                               function(VarUse){#(VarUse <- names(ClimNorm4RCPTmp[[1]])[1])
                                 app(do.call("c",
                                             lapply(ClimNorm4RCPTmp,
                                                    function(x){x[[VarUse]]})),mean)
                               }))
  names(ClimNormMn) <- names(ClimNorm4RCPTmp[[1]])
  # Estimate the SD of each band for the Climate Normal
  ClimNormSD <- do.call("c",
                        lapply(names(ClimNorm4RCPTmp[[1]]),
                               function(VarUse){#(VarUse <- names(ClimNorm4RCPTmp[[1]])[1])
                                 app(do.call("c",
                                             lapply(ClimNorm4RCPTmp,
                                                    function(x){x[[VarUse]]})),sd)
                               }))
  names(ClimNormSD) <- names(ClimNorm4RCPTmp[[1]])
  #Crop oceans and save the Climate Normal summary rasters
  ClimNormMn <- mask(ClimNormMn,BiomeBsLn)
  writeRaster(ClimNormMn,
              paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/ClimNormMn_1980-2010_",RCP,".tif"),
              overwrite = TRUE)
  ClimNormSD <- mask(ClimNormSD,BiomeBsLn)
  writeRaster(ClimNormSD,
              paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/ClimNormSD_1980-2010_",RCP,".tif"),
              overwrite = TRUE)
  rm(list=c("ClimNorm4RCPTmp","ClimNormMn","ClimNormSD"));gc()
}

