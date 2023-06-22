rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
require(tidyverse)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")

WGSRast <- rast(nrows=180, ncols=360) 
eck4Rast <- project(WGSRast,"+proj=eck4")
eck4Rast <- rast(extent=ext(eck4Rast),resolution=50000,crs="+proj=eck4")
rm(WGSRast);gc()

# Estimate Novelty based on Future to - CLimate Normal distance
#####
#####
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
  # Define the "BIOME" for each location
  if(!"WWF_BIOME_eck4_50km.tif"%in%dir("./Data/WWF-Biomes")){
    # Biomes Map
    BIOMES <- vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")
    BIOMES <- project(BIOMES,"+proj=eck4")
    BiomeBsLn <- rasterize(BIOMES,ClimNormMn[[1]],field = "BIOME")
    BiomeBsLn <- app(BiomeBsLn,function(x){ifelse(is.na(x),NA,ifelse(x>14,NA,x))})
    writeRaster(BiomeBsLn,
                "./Data/WWF-Biomes/WWF_BIOME_eck4_50km.tif")
  } else {
    BiomeBsLn <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_50km.tif")
  }
  
  #Crop oceans and save the Climate Normal summary rasters
  ClimNormMn <- mask(ClimNormMn,BiomeBsLn)
  writeRaster(ClimNormMn,
              paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM2/ClimNormMn_1980-2010_",RCP,".tif"),
              overwrite = TRUE)
  ClimNormSD <- mask(ClimNormSD,BiomeBsLn)
  writeRaster(ClimNormSD,
              paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM2/ClimNormSD_1980-2010_",RCP,".tif"),
              overwrite = TRUE)
  rm(list=c("ClimNorm4RCPTmp","ClimNormMn","ClimNormSD"));gc()
  #####
  #####
  # Estimate Novelty threshold based on Climate-Normal distance
  if(!paste0("AllModels_",RCP,"_TreshSumm.rds")%in%dir(paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM2/"))){
    
    # Estimate pairwise mahalanobis differences for the Climate-Normal period using a parellezed aprroach
    TimeCont <- Sys.time()
    sfInit(parallel=TRUE, cpus=200)
    sfExport("RCP")
    sfLibrary(terra)
    MDtreshDistList <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                                function(x,
                                         RCPUse = RCP){
                                  # Data frame with the Values
                                  ClimMn <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/ClimNormMn_1980-2010_",RCPUse,".tif"))
                                  ClimMn <- values(ClimMn, na.rm = T)
                                  # Data frame with the Covariance Matrix
                                  CoVarMtrx <- cov(ClimMn)
                                  # Estimate the mahalanobis Distance
                                  out <- mahalanobis(ClimMn,
                                                     ClimMn[x,],
                                                     CoVarMtrx)
                                  out <- data.frame(MnD = round(out,3))
                                  # Save the summary as a temp csv
                                  return(out)
                                })
    sfStop()
    # Merge the output in a single file
    MDtreshDist <- do.call("rbind",MDtreshDistList)
    # Estimate the mahalanobis distance based no-anlaogue threshold 
    MDtresh <- roc(object = MDtreshDist,
                   groups = values(BiomeBsLn, na.rm=T))
    rm(list = c("MDtreshDist","MDtreshDistList"));gc()
    Sys.time() - TimeCont
    # Estimate pairwise SED differences for the Climate-Normal period
    
    TimeCont <- Sys.time()
    sfInit(parallel=TRUE, cpus=200)
    sfExport("RCP")
    sfLibrary(terra)
    SEDtreshDistList <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                                 function(x,
                                          RCPUse = RCP){
                                   # Data frame with the Values
                                   ClimMn <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/ClimNormMn_1980-2010_",RCPUse,".tif"))
                                   ClimSD <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/ClimNormSD_1980-2010_",RCPUse,".tif"))
                                   TrgCellVals <- values(ClimMn,na.rm=T)[x,]
                                   # Estimate the Stdz Euc Distance Distance
                                   out <- as.numeric(values(sum(((ClimMn-TrgCellVals)^2)/ClimSD)^0.5,na.rm=T))
                                   out <- round(out,3)
                                   return(out)})
    sfStop()
    
    SEDtreshDist <- do.call("rbind",SEDtreshDistList)
    rm(SEDtreshDistList);gc()             
    # Estimate the SED distance based no-anlaogue threshold 
    SEDtresh <- roc(object = SEDtreshDist,
                    groups = values(BiomeBsLn, na.rm=T))
    rm(list = c("SEDtreshDist","SEDtreshDistList"))
    Sys.time() - TimeCont
    # final summary all Values
    Out.List <- list(RCP = RCP,
                     MDSummTresh = MDtresh$roc$Combined$optimal,
                     SEDSummTresh = SEDtresh$roc$Combined$optimal)
    #Save the Output
    saveRDS(Out.List,
            paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM2/AllModels_",RCP,"_TreshSumm.rds"))
    rm(list = c("Out.List","MDtresh","SEDtresh"));gc()
  }
}