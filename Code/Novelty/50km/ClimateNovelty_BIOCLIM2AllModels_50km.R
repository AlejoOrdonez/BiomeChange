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
              paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/ClimNormMn_1980-2010_",RCP,".tif"),
              overwrite = TRUE)
  ClimNormSD <- mask(ClimNormSD,BiomeBsLn)
  writeRaster(ClimNormSD,
              paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/ClimNormSD_1980-2010_",RCP,".tif"),
              overwrite = TRUE)
  rm(list=c("ClimNorm4RCPTmp","ClimNormMn","ClimNormSD"));gc()
  #####
  #####
  # Estimate Novelty threshold based on Climate-Normal distance
  if(!paste0("AllModels_",RCP,"_TreshSumm.rds")%in%dir(paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/"))){
    
    # Estimate pairwise mahalanobis differences for the Climate-Normal period using a parellezed aprroach
    sfInit(parallel=TRUE, cpus=100)
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
                                  return(out)
                                })
    sfStop()
    # Merge the list into a Table
    MDtreshDist <- do.call("cbind",MDtreshDistList)
    rm(MDtreshDistList);gc()
    # Estimate the mahalanobis distance based no-anlaogue threshold 
    MDtresh <- roc(object = MDtreshDist,
                   groups = values(BiomeBsLn, na.rm=T))
    rm(MDtreshDist);gc()
    
    # Estimate pairwise SED differences for the Climate-Normal period
    sfInit(parallel=TRUE, cpus=100)
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
                                   out <- data.frame(SED = round(out,3))
                                   return(out)})
    sfStop()
    # Merge the list into a Table
    SEDtreshDist <- do.call("cbind",SEDtreshDistList)
    rm(SEDtreshDistList);gc()             
    # Estimate the SED distance based no-anlaogue threshold 
    SEDtresh <- roc(object = SEDtreshDist,
                    groups = values(BiomeBsLn, na.rm=T))
    rm(SEDtreshDist);gc()
    
    # final summary all Values
    Out.List <- list(RCP = RCP,
                     MDSummTresh = MDtresh$roc$Combined$optimal,
                     SEDSummTresh = SEDtresh$roc$Combined$optimal)
    #Save the Output
    saveRDS(Out.List,
            paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/AllModels_",RCP,"_TreshSumm.rds"))
    rm(list = c("Out.List","MDtresh","SEDtresh"));gc()
  }
  
  #####
  #####
  # Estimate Novelty by comparing future climate to the Climate-Normal
  for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[1])
    if(!paste0("AllModels_",RCP,"_",YearUse,"_SEDminSumm.tif")%in%dir(paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/"))){
      # Load the Future points (and 11 yr period of conditions centred at the point of interest) data for an RCP
      RCPList <- lapply(c(YearUse-5):c(YearUse+5),
                        function(YearUseTmp){#(YearUse <- c(c(YearUse-9):YearUse)[1])
                          tmp <- rast(paste0("./Data/CMIP5/Processed/AllModelSumm1Arcmin/",RCP,"/BIOCLIM/AllModels_",RCP,"_",YearUse,".tif"),
                                      lyrs= BioclimVars)
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
                  paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/RCPFull_",YearUse,"_",RCP,".tif"),
                  overwrite = TRUE)
      rm(list = c("RCPFull","RCPList"));gc()
      # Estimate for each Future Clime ensemble, where is the closest analogue using the mahalanobis Distance
      
      sfInit(parallel=TRUE, cpus=100)
      sfExport("RCP")
      sfExport("YearUse")
      sfLibrary(terra)

      # Estimate the min mahalanobis distance per furture grid and save it as a temp CSV 
      MDmin <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                        function(x,
                                 RCPUse = RCP,
                                 YearFut = YearUse){
                          # Future condition
                          TrgCellVals <- values(rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/RCPFull_",YearFut,"_",RCPUse,".tif")),na.rm=T)[x,]
                          # Baseline Data frame with the Values
                          ClimMn <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/ClimNormMn_1980-2010_",RCPUse,".tif"))
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
                                            MD.min = round(min(MD.min,na.rm = TRUE),3), # Define the MDmin value of a future cell to all Climate normal cell
                                            Cell = sample(as.numeric(names(MD.min)[which(MD.min==min(MD.min,na.rm = TRUE))]),1))
                          # Save the summary as a temp CSV
                          data.table::fwrite(out,
                                             paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/Temp/Temp_MD_",x,".csv"),
                                             col.names = F)
                          gc()
                          return(TRUE)
                        })
      sfStop()
      
      # Merge the save CSVs into a single distance summary
      a <- paste0("/Users/au467796/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Results2/Novelty/AllModels_50km/",
                  RCP,"/BIOCLIM/")
      system(paste0("cat '",a,"Temp/Temp_MD_1.csv' '",
                    a,"Temp/Temp_MD_2.csv' >'",a,"All_MDmin.csv'"))
      Brq <- round(seq(2,
                       dim(values(BiomeBsLn,na.rm=T))[1],
                       length.out=25))
      for(j in 1:24){#j<-24
        i <- c(c(Brq[j]+1):c(Brq[j+1]))
        system(paste0 ("cat ",
                       paste0("'",a,"Temp/Temp_MD_",i,".csv'",collapse = " "),
                       " >> '",a,"All_MDmin.csv'"))
        }
      rm(list=c("j","i","Brq"));gc()

      # Load the Appended min mahalanobis distance table
      MDmin <- read.csv(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/All_MDmin.csv"),
                     header = F,
                     col.names = c("ID","MD.min","Cell"))
      
      # Turn the min mahalanobis distance table into a SpatRast
      MDminrast <- rast(BiomeBsLn, nlyrs=2)
      MDminrast[!is.na(BiomeBsLn[])] <- MDmin[,-1] #do.call("rbind",MDmin)
      system(paste0("rm -r '",a,"Temp/'*"))
      
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
                  paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/AllModels_",RCP,"_",YearUse,"_MDminSumm.tif"),
                  overwrite = TRUE)
      rm(list=c("a","MDmin","MDMinDistinKm","MDminSumm"));gc()
      
      
      # Estimate for each Future Clime ensemble, where is the closest analogue using the Standarized Euclidean Distance
      TimeIn2 <- Sys.time()
      sfInit(parallel=TRUE, cpus=100)
      sfExport("RCP")
      sfExport("YearUse")
      sfLibrary(terra)
      SEDMin <- sfLapply(1:dim(values(BiomeBsLn,na.rm=T))[1],
                         function(x,
                                  RCPUse = RCP,
                                  YearFut = YearUse){
                           # Future condition
                           TrgCellVals <- values(rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/RCPFull_",YearFut,"_",RCPUse,".tif")),na.rm=T)[x,]
                           # Baseline Data frame with the Values
                           ClimMn <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/ClimNormMn_1980-2010_",RCPUse,".tif"))
                           ClimSD <- rast(paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/ClimNormSD_1980-2010_",RCP,".tif"))
                           # Estimate the Standarzied euclidean Distance
                           SEDRast <- sum(((ClimMn-TrgCellVals)^2)/ClimSD)^0.5
                           # compile the Min distance and Cell ID
                           out <- data.frame(ID = x,
                                             SED.Min = round(minmax(SEDRast)['min',],3), # Define the MDmin value of a future cell to all Climate normal cell
                                             Cell = sample(which(values(SEDRast)==(minmax(SEDRast)['min',])),1))
                           # Save the summary as a temp CSV
                           data.table::fwrite(out,
                                              paste0("./Results2/Novelty/AllModels_50km/",RCPUse,"/BIOCLIM/Temp/Temp_SED_",x,".csv"),
                                              col.names = F)
                           gc()
                           return(TRUE)
                         })
      sfStop()
      # Merge the save CSVs into a single distance summary
      a <- paste0("/Users/au467796/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Results2/Novelty/AllModels_50km/",
                  RCP,"/BIOCLIM/")
      system(paste0("cat '",a,"Temp/Temp_SED_1.csv' '",
                    a,"Temp/Temp_SED_2.csv' >'",a,"All_SEDmin.csv'"))
      Brq <- round(seq(2,
                       dim(values(BiomeBsLn,na.rm=T))[1],
                       length.out=25))
      for(j in 1:24){#j<-24
        i <- c(c(Brq[j]+1):c(Brq[j+1]))
        system(paste0 ("cat ",
                       paste0("'",a,"Temp/Temp_SED_",i,".csv'",collapse = " "),
                       " >> '",a,"All_SEDmin.csv'"))
      }
      rm(list=c("j","i","Brq"));gc()
      
      # Load the Appended min SED distance table
      SEDMin <- read.csv(paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/All_SEDmin.csv"),
                        header = F,
                        col.names = c("ID","SED.Min","Cell"))
      # Turn the min SED distance table into a SpatRast
      SEDminrast <- rast(BiomeBsLn, nlyrs=2)
      SEDminrast[!is.na(BiomeBsLn[])] <- SEDMin[,-1]
      names(SEDminrast) <- c("SED.Min","Cell")
      #system(paste0("rm -r '",a,"Temp/'*"))
      #file.remove(dir(paste0(a,"Temp"),
      #                full.names = T))
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
                  paste0("./Results2/Novelty/AllModels_50km/",RCP,"/BIOCLIM/AllModels_",RCP,"_",YearUse,"_SEDminSumm.tif"),
                  overwrite = TRUE)      
      rm(list=c("a","SEDmin","SEDMinDistinKm","SEDminSumm"));gc()
      Sys.time()-TimeIn2
    }
  }
  gc()
}

