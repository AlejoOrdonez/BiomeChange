rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
require(tidyverse)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")

# Create a 100 x100 km raster
WGSRast <- rast(nrows=180, ncols=360) 
eck4Rast <- project(WGSRast,"+proj=eck4")
eck4Rast <- rast(extent=ext(eck4Rast),resolution=50000,crs="+proj=eck4")
rm(WGSRast);gc()

# Load the BIOME
BiomeBsLn <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_50km.tif")

# Estimate Novelty threshold based on Climate-Normal distance
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[4])
  #####
  if(!paste0("AllModels_",RCP,"_TreshSumm.rds")%in%dir(paste0("./Results2/Novelty/AllModels_50km/",RCP,"/koeppen_geiger/"))){
    # Estimate pairwise mahalanobis differences for the Climate-Normal period using a parellezed aprroach
    TimeCont2 <- Sys.time()
    sfInit(parallel=TRUE, cpus=100)
    sfExport("RCP")
    sfLibrary(terra)
    MDtreshDistList <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                                function(x,
                                         RCPUse = RCP){
                                  # Data frame with the Values
                                  ClimMn <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/koeppen_geiger/ClimNormMn_1980-2010_",RCPUse,".tif"))
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
    MDtreshDist <- do.call("cbind",MDtreshDistList)
    rm(MDtreshDistList);gc()
    # Estimate the mahalanobis distance based no-anlaogue threshold 
    MDtresh <- roc(object = MDtreshDist,
                   groups = values(BiomeBsLn, na.rm=T))
    rm(list = c("MDtreshDist"));gc()
    # Estimate pairwise SED differences for the Climate-Normal period
    sfInit(parallel=TRUE, cpus=100)
    sfExport("RCP")
    sfLibrary(terra)
    SEDtreshDistList <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                                 function(x,
                                          RCPUse = RCP){
                                   # Data frame with the Values
                                   ClimMn <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/koeppen_geiger/ClimNormMn_1980-2010_",RCPUse,".tif"))
                                   ClimSD <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/koeppen_geiger/ClimNormSD_1980-2010_",RCPUse,".tif"))
                                   TrgCellVals <- values(ClimMn,na.rm=T)[x,]
                                   # Estimate the Stdz Euc Distance Distance
                                   out <- as.numeric(values(sum(((ClimMn-TrgCellVals)^2)/ClimSD)^0.5,na.rm=T))
                                   if(length(out)!=54200){
                                     out <- as.numeric(values(sum(((ClimMn-TrgCellVals)^2)/ClimSD)^0.5))
                                     out <- out[!is.na(BiomeBsLn[])]
                                     out[is.na(out)] <- min(out,na.rm=T)
                                   }
                                   out <- data.frame(SED = round(out,3))
                                   return(out)})
    sfStop()
    SEDtreshDist <- do.call("cbind",SEDtreshDistList)
    rm(SEDtreshDistList);gc()             
    # Estimate the SED distance based no-anlaogue threshold 
    SEDtresh <- roc(object = SEDtreshDist,
                    groups = values(BiomeBsLn, na.rm=T))
    rm(list = c("SEDtreshDist"))
    # final summary all Values
    Out.List <- list(RCP = RCP,
                     MDSummTresh = MDtresh$roc$Combined$optimal,
                     SEDSummTresh = SEDtresh$roc$Combined$optimal)
    #Save the Output
    saveRDS(Out.List,
            paste0("./Results2/Novelty/AllModels_50km/",RCP,"/koeppen_geiger/AllModels_",RCP,"_TreshSumm.rds"))
    rm(list = c("Out.List","MDtresh","SEDtresh"));gc()
    print(Sys.time() - TimeCont2)
    
  }
}