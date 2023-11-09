rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
require(tidyverse)

#####
# Estimate Novelty based on Future to - CLimate Normal distance
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[1])
  for(ModelUse in c("bcc-csm1-1", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0",
                    "GISS-E2-H", "GISS-E2-R", "HadGEM2-ES","MPI-ESM-LR", "NorESM1-M")){#(ModelUse <- "NorESM1-M")
    # Load the Climate Normal Data
    ClimMn_50x50 <- c(rast(dir("./Data/Raw/By_Model/Historical",pattern=ModelUse, full.names = T)[1]),
                      rast(dir("./Data/Raw/By_Model/Historical",pattern=ModelUse, full.names = T)[2]))
    # Estimate Novelty by comparing future climate to the Climate-Normal
    for(YearUse in seq(2100,2300,by=50)){#(YearUse <- seq(2100,2300,by=50)[3])
      # Ensure there is future data
      if(length(dir(paste0("./",RCP),pattern = paste0(ModelUse,"_",YearUse-15)))!=0){
        # Ensure there is future data and that the Disimilarity has not been estimated
        TimeIn2 <- Sys.time()
        # Load the Future CLimatic conditions 
        RCPFull_50x50 <- c(rast(dir(paste0("./Data/Raw/By_Model/",RCP),pattern = paste0(ModelUse,"_",YearUse-15),full.names = T)[1]),
                           rast(dir(paste0("./Data/Raw/By_Model/",RCP),pattern = paste0(ModelUse,"_",YearUse-15),full.names = T)[2]))
        # Build a values Data frame for the climate normal [1:8] and the Future conditions [9:16]
        Normal_Future_Climate_rast <- values(c(ClimMn_50x50,
                                               RCPFull_50x50),
                                             na.rm=T)
        # Scale the climate values Data frame
        Normal_Future_Climate_rastScl <- scale(Normal_Future_Climate_rast)
        
        # Define Cells with data
        cell_with_data <- which(complete.cases(values(c(ClimMn_50x50[[1]],RCPFull_50x50[[1]]))))
        
        # Estimate for each Future Clime ensemble, the distance (in KM)  to the closest analogue using the mahalanobis Distance
        CordsNormal <- crds(ClimMn_50x50, df = F) # get the coordinates of a cell
        
        # Estimate for each Future Clime ensemble, where is the closest analogue using the mahalanobis Distance
        sfInit(parallel=TRUE, cpus=20)
        sfExport("Normal_Future_Climate_rastScl")
        sfExport("CordsNormal")
        # Estimate the min mahalanobis distance per furture grid and save it as a temp CSV 
        MDMin_list <- sfLapply(1:dim(Normal_Future_Climate_rastScl)[1],
                               function(x){#(x<-30947)
                                 # Future condition
                                 TrgCellVals <- Normal_Future_Climate_rastScl[x,9:16]
                                 # Baseline Data frame with the Values
                                 ClimMnTbl <- Normal_Future_Climate_rastScl[,1:8]
                                 # Data frame with the Covariance Matrix
                                 CoVarMtrx <- cov(ClimMnTbl)
                                 # Estimate the mahalanobis Distance
                                 MD.min <- mahalanobis(ClimMnTbl,
                                                       TrgCellVals,
                                                       CoVarMtrx)
                                 # Estimate for each Future Clime ensemble, the distance (in KM)  to the closest analogue using the mahalanobis Distance
                                 Dist_to_MD.min <- as.vector(terra::distance(x = matrix(CordsNormal[x,],ncol=2),
                                                                             y = matrix(CordsNormal[which(MD.min==min(MD.min,na.rm = TRUE)),],ncol=2),
                                                                             lonlat = F)/1000)
                                 # compile the Min distance and Cell ID
                                 out <- data.frame(ID = x,
                                                   MD.min = round(min(MD.min,na.rm = TRUE),3), # Define the MDMin value of a future cell to all Climate normal cell
                                                   Cell = which(MD.min==min(MD.min,na.rm = TRUE))[which(Dist_to_MD.min==min(Dist_to_MD.min))][1], # Define the Cell with the closest analogue
                                                   DistinKm = min(Dist_to_MD.min)  # Define the distance to the Closest analogue
                                 )
                                 return(out)
                               })
        sfStop()
        # Load the Appended min mahalanobis distance table
        MDMin <- do.call("rbind",MDMin_list)
        rm(list = c("MDMin_list","Normal_Future_Climate_rastScl","Normal_Future_Climate_rast","RCPFull_50x50","RCPFullPrj","RCPFull"));gc()
        # Turn the min mahalanobis distance, Cell ID and distance table into a SpatRast
        MDminSumm <- rast(CordsNormal[[1]], nlyrs=3)
        MDminSumm[cell_with_data] <- MDMin[,-1] #do.call("rbind",MDMin)
        names(MDminSumm) <- c("MD.min", "Cell", "DistinKm")
        writeRaster(MDminSumm,
                    paste0("./Results/Novelty/By_Model/",RCP,"/MDminSumm_",ModelUse,"_",RCP,"_",YearUse,".tif"),
                    overwrite = TRUE)
        rm(list=c("MDminSumm", "CordsNormal","cell_with_data"));gc()
        print(Sys.time()-TimeIn2)
      }
      print(YearUse)
    }
    print(ModelUse)
    rm(list=c("ClimMn_50x50"));gc()
  }
  print(RCP)  
}
