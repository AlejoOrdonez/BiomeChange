rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
setwd("~/Documents/GitHub/BiomeChange")


# Define the "BIOME" for each location
BiomeBsLn <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_50km.tif")

# Estimate the baseline mean and SD of climate variables
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[1])
  # Estimate Novelty by comparing future climate to the Climate-Normal
  for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[2])
    if(!paste0("AllModels_",RCP,"_",YearUse,"_SEDminSumm.tif")%in%dir(paste0("./Results/Novelty/Mean_All_Models/",RCP))){
      # Estimate for each Future Clime ensemble, where is the closest analogue using the mahalanobis Distance
      TimeIn <- Sys.time()
      sfInit(parallel=TRUE, cpus=6)
      sfExport("RCP")
      sfExport("YearUse")
      sfLibrary(terra)
      MDmin <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                        function(x,
                                 RCPUse = RCP,
                                 YearFut = YearUse){
                          # Load the Future points (and 11 yr period of conditions centred at the point of interest) data for an RCP
                          TrgCellVals <- values(rast(paste0("./Data/Raw/Mean_All_Models/",RCP,"/RCPFull_",YearUse,"_",RCP,".tif")),na.rm=T)[x,]
                          # Load the 1980 to 2010 Historical data for an RCP to Create a climate normal raster for each evaluated variable  
                          ClimNormMn <- rast(paste0("./Data/Raw/Mean_All_Models/",RCP,"/ClimNormMn_1980-2010_",RCP,".tif"))
                          ClimMnTbl <- values(ClimNormMn,na.rm=T)
                          rownames(ClimMnTbl) <- which(!is.na(values(ClimNormMn[[1]])))
                          # Data frame with the Covariance Matrix
                          CoVarMtrx <- cov(na.omit(ClimMnTbl))
                          # Estimate the mahalanobis Distance
                          MD.min <- mahalanobis(ClimMnTbl,
                                                TrgCellVals,
                                                CoVarMtrx)
                          out <- c(MD.min = round(min(MD.min,na.rm = TRUE),3), # Define the MDmin value of a future cell to all Climate normal cell
                                   Cell = as.numeric(names(MD.min)[which(MD.min==min(MD.min,na.rm = TRUE))]))
                          return(out)
                        })
      sfStop()
      MDminrast <- rast(BiomeBsLn, nlyrs=2)
      MDminrast[!is.na(BiomeBsLn[])] <- do.call("rbind",MDmin)
      names(MDminrast) <- c("MD.min","Cell")
      Sys.time() - TimeIn
      
      # Estimate for each Future Clime ensemble, the distance (in KM)  to the closest analogue using the mahalanobis Distance
      TimeIn <- Sys.time()
      CordsAll <- crds(BiomeBsLn, df = TRUE,na.rm = FALSE) # get the coordinates
      MDMinDistinKm <- rast(BiomeBsLn, nlyrs = 3) # Make am empty SpatRaster File to summarize values
      values(MDMinDistinKm) <- cbind(CordsAll,
                                     values(MDminrast[[2]])) # Add values (Coordinates) and the closest Cell
      names(MDMinDistinKm) <- c("x","y","ID") # Renames so the distance function works (positions need to be x/y OR lat/lon)
      # Estimate the distance of each cell to the Closest analogue
      MDMinDistinKm <- app(MDMinDistinKm,
                           function(x){#(x<-values(MDMinDistinKm)[1000,])
                             if(!is.na(x[3])){
                               Dist <- as.numeric(min(terra::distance(x = as.matrix(t(x[1:2])),
                                                                      y = as.matrix(CordsAll[x[3],]),
                                                                      lonlat = F)))/1000  
                             } else{
                               Dist <- NA
                             }
                             return(Dist)
                           })
      Sys.time() - TimeIn
      # Make a summary for mahalanobis distance estimates  
      MDminSumm <- c(MDminrast,MDMinDistinKm)
      names(MDminSumm) <- c(names(MDminrast),"DistinKm")
      writeRaster(MDminSumm,
                  paste0("./Results/Novelty/Mean_All_Models/",RCP,"/AllModels_",RCP,"_",YearUse,"_MDminSumm.tif"),
                  overwrite = TRUE)
      rm(list=c("MDmin","MDMinDistinKm","MDminSumm"));gc()
    }
  }

}
