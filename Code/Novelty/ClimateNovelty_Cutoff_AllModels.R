rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
setwd("~/Documents/GitHub/BiomeChange")

# Load the BIOME
BiomeBsLn <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_50km.tif")

# Estimate Novelty threshold based on Climate-Normal distance
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[1])
  #####
  if(!paste0("AllModels_",RCP,"_TreshSumm.rds")%in%dir(paste0("./Results/Novelty/Mean_All_Models/",RCP))){
    # Estimate pairwise mahalanobis differences for the Climate-Normal period using a parellezed aprroach
    TimeCont2 <- Sys.time()
    sfInit(parallel=TRUE, cpus=100)
    sfExport("RCP")
    sfLibrary(terra)
    MDtreshDistList <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                                function(x,
                                         RCPUse = RCP){
                                  # Load the 1980 to 2010 Historical data for an RCP to Create a climate normal raster for each evaluated variable  
                                  ClimNormMn <- rast(paste0("./Data/Raw/Mean_All_Models/",RCP,"/ClimNormMn_1980-2010_",RCP,".tif"))
                                  ClimMn <- values(ClimNormMn, na.rm = T)
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
    # final summary all Values
    Out.List <- list(RCP = RCP,
                     MDSummTresh = MDtresh$roc$Combined$optimal)
    #Save the Output
    saveRDS(Out.List,
            paste0("./Results/Novelty/Mean_All_Models/",RCP,"/AllModels_",RCP,"_TreshSumm.rds"))
    rm(list = c("Out.List","MDtresh"));gc()
    print(Sys.time() - TimeCont2)
  }
}