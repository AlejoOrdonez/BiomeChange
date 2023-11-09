rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
require(tidyverse)

# Load the BIOME
BiomeBsLn <- rast("~/Documents/GitHub/BiomeChange/Data/WWF-Biomes/WWF_BIOME_eck4_50km.tif")

# Estimate Novelty threshold based on Climate-Normal distance
for(ModelUse in c("bcc-csm1-1", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0",
                    "GISS-E2-H", "GISS-E2-R", "HadGEM2-ES","MPI-ESM-LR", "NorESM1-M")){#(ModelUse <- "bcc-csm1-1")
  # Estimate pairwise mahalanobis differences for the Climate-Normal period using a parellezed aprroach
  # Load the Climate Normal Data
  ClimMn <- c(rast(dir("~/Documents/GitHub/BiomeChange/Data/Raw/Historical",
                       pattern=ModelUse, full.names = T)[1]),
              rast(dir("~/Documents/GitHub/BiomeChange/Data/Raw/Historical",
                       pattern=ModelUse, full.names = T)[2]))

  # Scale the climate mean
  ClimMnScle <- scale(values(ClimMn_50x50, na.rm = T))
  # Data frame with the Covariance Matrix
  CoVarMtrx <- cov(ClimMnScle)
  # estimate all pairwise distances
  TimeCont1 <- Sys.time()
  sfInit(parallel=TRUE, cpus=20)
  sfExport("ClimMnScle")
  sfExport("CoVarMtrx")
  sfLibrary(terra)
  MDtreshDistList <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                              function(x){#(x<-1)
                                # Estimate the mahalanobis Distance
                                out <- mahalanobis(ClimMnScle,
                                                   ClimMnScle[x,],
                                                   CoVarMtrx)
                                out <- data.frame(MnD = round(out,3))
                                # Save the summary as a temp csv
                                return(out)
                              })
  sfStop()
  Sys.time()-TimeCont1
  # Merge the output in a single file
  MDtreshDist <- do.call("cbind",MDtreshDistList)
  rm(MDtreshDistList);gc()
  # Estimate the mahalanobis distance based no-anlaogue threshold 
  MDtresh <- roc(object = MDtreshDist,
                 groups = values(BiomeBsLn, na.rm=T))
  rm(list = c("MDtreshDist"));gc()
  # final summary all Values
  Out.List <- list(Model = ModelUse,
                   MDSummTresh = MDtresh$roc$Combined$optimal)
  print(Sys.time() - TimeCont2)
  Out.List
  #Save the Output
  saveRDS(Out.List,
          paste0("./Results/Novelty/Tresholds/",
                 "TreshSumm_",ModelUse,".rds"))
  rm(list = c("Out.List","MDtresh","SEDtresh"));gc()
  print(Sys.time() - TimeCont2)
}
  
  
