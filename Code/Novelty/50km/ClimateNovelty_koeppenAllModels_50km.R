rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")

# Create a 100 x100 km raster
WGSRast <- rast(nrows=180, ncols=360) 
eck4Rast <- project(WGSRast,"+proj=eck4")
eck4Rast <- rast(extent=ext(eck4Rast),resolution=50000,crs="+proj=eck4")
rm(WGSRast);gc()

# Load the bimome raster
BiomeBsLn <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_50km.tif")

# Estimate Novelty based on Future to - Climate Normal distance
#####
#####
# Estimate the baseline mean and SD of climate variables
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[2])
  #####
  # Estimate Novelty by comparing future climate to the Climate-Normal
  for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[5])
    TimeIn2 <- Sys.time()
    # Load the Future points (and 11 yr period of conditions centred at the point of interest) data for an RCP
    RCPList <- lapply(c(YearUse-5):c(YearUse+5),
                      function(YearUseTmp){#(YearUse <- c(c(YearUse-9):YearUse)[1])
                        tmp <- rast(paste0("./Data/CMIP5/Processed/AllModelSumm1Arcmin/",RCP,"/koeppen_geiger/AllModels_",RCP,"_",YearUse,".tif"))
                        tmp <- project(tmp,"+proj=eck4")
                        resample(tmp,
                                 eck4Rast,
                                 method = "near")
                      })
    # Estimate the mean of each band for the 
    RCPFull <- Reduce("+",RCPList)/length(RCPList)
    #Crop oceans
    RCPFull <- mask(RCPFull,BiomeBsLn)
    writeRaster(RCPFull,
                paste0("./Results2/Novelty/AllModels_50km/",RCP,"/koeppen_geiger/RCPFull_",YearUse,"_",RCP,".tif"),
                overwrite = TRUE)
    rm(list = c("RCPFull","RCPList"));gc()
    # Estimate for each Future Clime ensemble, where is the closest analogue using the mahalanobis Distance
    if(!paste0("AllModels_",RCP,"_",YearUse,"_MDminSumm.tif")%in%dir(paste0("./Results2/Novelty/AllModels_50km/",RCP,"/koeppen_geiger/"))){
      sfInit(parallel=TRUE, cpus=100)
      sfExport("RCP")
      sfExport("YearUse")
      sfLibrary(terra)
      
      # Estimate the min mahalanobis distance per furture grid and save it as a temp CSV 
      MDMin <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                        function(x,
                                 RCPUse = RCP,
                                 YearFut = YearUse){
                          # Future condition
                          TrgCellVals <- values(rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/koeppen_geiger/RCPFull_",YearFut,"_",RCPUse,".tif")),na.rm=T)[x,]
                          # Baseline Data frame with the Values
                          ClimMn <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/koeppen_geiger/ClimNormMn_1980-2010_",RCPUse,".tif"))
                          ClimMnTbl <- values(ClimMn,na.rm=T)
                          rownames(ClimMnTbl) <- which(!is.na(values(ClimMn[[1]])))
                          # Data frame with the Covariance Matrix
                          CoVarMtrx <- cov(na.omit(ClimMnTbl))
                          # Estimate the mahalanobis Distance
                          MD.min <- mahalanobis(ClimMnTbl,
                                                TrgCellVals,
                                                CoVarMtrx)
                          # compile the Min distance and Cell ID
                          out <- data.frame(ID = x,
                                            MD.min = round(min(MD.min,na.rm = TRUE),3), # Define the MDMin value of a future cell to all Climate normal cell
                                            Cell = sample(as.numeric(names(MD.min)[which(MD.min==min(MD.min,na.rm = TRUE))]),1))
                          return(out)
                        })
      sfStop()
      # Load the Appended min mahalanobis distance table
      MDMin <- do.call("rbind",MDMin)
      
      # Turn the min mahalanobis distance table into a SpatRast
      MDminrast <- rast(BiomeBsLn, nlyrs=2)
      MDminrast[!is.na(BiomeBsLn[])] <- MDMin[,-1] #do.call("rbind",MDMin)
      
      # Estimate for each Future Clime ensemble, the distance (in KM)  to the closest analogue using the mahalanobis Distance
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
      
      # Make a summary for mahalanobis distance estimates  
      MDminSumm <- c(MDminrast,MDMinDistinKm)
      names(MDminSumm) <- c(names(MDminrast),"DistinKm")
      writeRaster(MDminSumm,
                  paste0("./Results2/Novelty/AllModels_50km/",RCP,"/koeppen_geiger/AllModels_",RCP,"_",YearUse,"_MDminSumm.tif"),
                  overwrite = TRUE)
      rm(list=c("MDMin","MDMinDistinKm","MDminSumm","MDminrast","CordsAll"));gc()
    }      
    if(!paste0("AllModels_",RCP,"_",YearUse,"_SEDminSumm.tif")%in%dir(paste0("./Results2/Novelty/AllModels_50km/",RCP,"/koeppen_geiger/"))){
      # Estimate for each Future Clime ensemble, where is the closest analogue using the Standarized Euclidean Distance
      
      sfInit(parallel=TRUE, cpus=100)
      sfExport("RCP")
      sfExport("YearUse")
      sfLibrary(terra)
      SEDMin <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                         function(x,
                                  RCPUse = RCP,
                                  YearFut = YearUse){
                           # Future condition
                           TrgCellVals <- values(rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/koeppen_geiger/RCPFull_",YearFut,"_",RCPUse,".tif")),na.rm=T)[x,]
                           # Baseline Data frame with the Values
                           ClimMn <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/koeppen_geiger/ClimNormMn_1980-2010_",RCPUse,".tif"))
                           ClimSD <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCP,"/koeppen_geiger/ClimNormSD_1980-2010_",RCP,".tif"))
                           # Estimate the Standarzied euclidean Distance
                           SEDRast <- sum(((ClimMn-TrgCellVals)^2)/ClimSD)^0.5
                           # compile the Min distance and Cell ID
                           out <- data.frame(ID = x,
                                             SED.Min = round(minmax(SEDRast)['min',],3), # Define the MDMin value of a future cell to all Climate normal cell
                                             Cell = sample(which(values(SEDRast)==(minmax(SEDRast)['min',])),1))
                           return(out)
                         })
      sfStop()
      SEDMin <- do.call("rbind",SEDMin)
      # Turn the min SED distance table into a SpatRast
      SEDminrast <- rast(BiomeBsLn, nlyrs=2)
      SEDminrast[!is.na(BiomeBsLn[])] <- SEDMin[,-1]
      names(SEDminrast) <- c("SED.Min","Cell")
      # Estimate for each Future Clime ensemble, the distance (in KM)  to the closest analogue using the Standarized Euclidean Distance
      CordsAll <- crds(BiomeBsLn, df = TRUE,na.rm = FALSE) # get the coordinates
      SEDMinDistinKm <- rast(BiomeBsLn, nlyrs = 3)  # Make am empty SpatRaster File to summarize values
      values(SEDMinDistinKm) <- cbind(CordsAll,values(SEDminrast[[2]])) # Add values (Coordinates) and the closest Cell
      names(SEDMinDistinKm) <- c("x","y","ID") # Renames so the distance function works (positions need to be x/y OR lat/lon)
      # Estimate the distance of each cell to the Closest analogue
      SEDMinDistinKm<- app(SEDMinDistinKm,
                           function(x){#(x<-values(SEDMinDistinKm)[1000,])
                             if(!is.na(x[3])){
                               Dist <- as.numeric(min(terra::distance(x = as.matrix(t(x[1:2])),
                                                                      y = as.matrix(CordsAll[x[3],]),
                                                                      lonlat = F)))/1000  
                             } else{
                               Dist <- NA
                             }
                             return(Dist)
                           })
      
      # Make a summary for SED distance estimates  
      SEDminSumm <- c(SEDminrast,SEDMinDistinKm)
      names(SEDminSumm) <- c(names(SEDminrast),"DistinKm")
      writeRaster(SEDminSumm,
                  paste0("./Results2/Novelty/AllModels_50km/",RCP,"/koeppen_geiger/AllModels_",RCP,"_",YearUse,"_SEDminSumm.tif"),
                  overwrite = TRUE)      
      rm(list=c("SEDMin","SEDMinDistinKm","SEDminSumm","SEDminrast","CordsAll"));gc()
    }
    print(Sys.time()-TimeIn2)
    print(YearUse)
  }
  gc()
}

