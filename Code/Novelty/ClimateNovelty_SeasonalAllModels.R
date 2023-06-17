rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")

# Biomes Map
BIOMES <- vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")

# Estimate Novelty based on Future to - CLimate Normal distance
#####
#####
  # Estimate the baseline mean and SD of climate variables
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[2])
  # Load the 1980 to 2010 Historical data for an RCP to Create a climate normal raster for each evaluated variable  
  ClimNorm4RCPTmp <- lapply(1980:2010,
                            function(YearUse){#(YearUse <- c(1980:2010)[1])
                              c(rast(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/Seasonal/AllModels_",RCP,"_",YearUse,"_P.Seasonal.tif")),
                                rast(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/Seasonal/AllModels_",RCP,"_",YearUse,"_Tm.Seasonal.tif")))
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
  #Crop oceans
  ClimNormMn <- mask(ClimNormMn,BIOMES[BIOMES$BIOME<15,])
  
  # Estimate the SD of each band for the Climate Normal
  ClimNormSD <- do.call("c",
                        lapply(names(ClimNorm4RCPTmp[[1]]),
                               function(VarUse){#(VarUse <- names(ClimNorm4RCPTmp[[1]])[1])
                                 app(do.call("c",
                                             lapply(ClimNorm4RCPTmp,
                                                    function(x){x[[VarUse]]})),sd)
                               }))
  names(ClimNormSD) <- names(ClimNorm4RCPTmp[[1]])
  #Crop oceans
  ClimNormSD <- mask(ClimNormSD,BIOMES[BIOMES$BIOME<15,])
#####
#####
  # Estimate Novelty threshold based on Climate-Normal distance
  if(!paste0("AllModels_",RCP,"_TreshSumm.rds")%in%dir(paste0("./Results/Novelty/AllModels/",RCP,"/Seasonal/"))){
    # Define the "BIOME" for each location
    BiomeBsLn <- terra::extract(BIOMES[,"BIOME"],
                                crds(ClimNormMn[[1]],df=T))
    TimeIn <- Sys.time()
    ClimNormSD[BiomeBsLn$id.y[BiomeBsLn$BIOME>14]]<-NA
    ClimNormMn[BiomeBsLn$id.y[BiomeBsLn$BIOME>14]]<-NA
    BiomeBsLn <- BiomeBsLn[BiomeBsLn$BIOME<15,]
    
    # Estimate pairwise mahalanobis differences for the Climate-Normal period
    MDtreshDist <- apply(values(ClimNormMn, na.rm = T),
                         1,
                         function(TrgCellVals){#(TrgCellVals<-values(ClimNormMn, na.rm = T)[1,])
                           mahalanobis(values(ClimNormMn, na.rm=T),
                                       TrgCellVals,
                                       cov(values(ClimNormMn, na.rm=T)))
                         })
    # Estimate the mahalanobis distance based no-anlaogue threshold 
    MDtresh <- roc(object = MDtreshDist,
                   groups = BiomeBsLn$BIOME)
    rm(MDtreshDist);gc()
    Sys.time()-TimeIn
    # Estimate pairwise SED differences for the Climate-Normal period
    TimeIn <- Sys.time()
    SEDtreshDist <- apply(values(ClimNormMn),
                          1,
                          function(TrgCellVals){#(TrgCellVals<-values(ClimNormMn)[1,])
                            as.numeric(values(sum(((ClimNormMn-TrgCellVals)^2)/ClimNormSD)^0.5))
                          })
    # Estimate the SED distance based no-anlaogue threshold 
    SEDtresh <- roc(object = SEDtreshDist,
                    groups = BiomeBsLn$BIOME)
    rm(SEDtreshDist);gc()
    Sys.time()-TimeIn
    # final summary all Values
    Out.List <- list(RCP = RCP,
                     MDSummTresh = MDtresh$roc$Combined$optimal,
                     SEDSummTresh = SEDtresh$roc$Combined$optimal)
    #Save the Output
    saveRDS(Out.List,
            paste0("./Results/Novelty/AllModels/",RCP,"/Seasonal/AllModels_",RCP,"_TreshSumm.rds"))
  }
  rm(list = ls()[!ls()%in%c("RCP","BIOMES",
                            "MDtresh","SEDtresh","Out.List",
                            "ClimNormMn","ClimNormSD")]);gc()
#####
#####
  # Estimate Novelty by comparing future climate to the Climate-Normal
  for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[2])
    if(!paste0("AllModels_",RCP,"_",YearUse,"_SEDminSumm.tif")%in%dir(paste0("./Results/Novelty/AllModels/",RCP,"/Seasonal/"))){
      # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
      RCPList <- lapply(c(YearUse-9):YearUse,
                        function(YearUseTmp){#(YearUseTmp <- c(c(YearUse-9):YearUse)[1])
                          c(rast(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/Seasonal/AllModels_",RCP,"_",YearUseTmp,"_P.Seasonal.tif")),
                            rast(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/Seasonal/AllModels_",RCP,"_",YearUseTmp,"_Tm.Seasonal.tif")))
                          
                        })
      # Estimate the mean of each band for the 
      RCPFull <- Reduce("+",RCPList)/length(RCPList)
      #Crop oceans
      RCPFull <- mask(RCPFull,BIOMES[BIOMES$BIOME<15,])
      # Estimate for each Future Clime ensemble, where is the closest analogue using the mahalanobis Distance
      MDmin <- app(RCPFull,
                   function(TrgCellVals){#(TrgCellVals<-RCPFull[][1,])
                     if(!is.na(TrgCellVals[1])){
                       MD.min <- mahalanobis(values(ClimNormMn),
                                             TrgCellVals,
                                             cov(values(ClimNormMn, na.rm=T)))
                       out <- c(MD.min=min(MD.min,na.rm = TRUE), # Define the MDmin value of a future cell to all Climate normal cell
                                Cell=which(MD.min==min(MD.min,na.rm = TRUE)))
                     } else{
                       out <- c(MD.min=NA, # Define the MDmin value of a future cell to all Climate normal cell
                                Cell=NA)
                     }
                     return(out)
                   })
      # Estimate for each Future Clime ensemble, the distance (in KM)  to the closest analogue using the mahalanobis Distance
      CordsAll <- crds(RCPFull, df = TRUE,na.rm = FALSE) # get the coordinates
      MDMinDistinKm <- rast(ClimNormMn[[1:3]]) # Make am empty SpatRaster File to summarize values
      values(MDMinDistinKm) <- cbind(CordsAll,values(MDmin[[2]])) # Add values (Coordinates) and the closest Cell
      names(MDMinDistinKm) <- c("x","y","ID") # Renames so the distance function works (positions need to be x/y OR lat/lon)
      # Estimate the distance of each cell to the Closest analogue
      MDMinDistinKm <- app(MDMinDistinKm,
                           function(x){#(x<-values(MDMinDistinKm)[1000,])
                             if(!is.na(x[3])){
                               Dist <- as.numeric(min(terra::distance(x = as.matrix(t(x[1:2])),
                                                                      y = as.matrix(CordsAll[x[3],]),
                                                                      lonlat = T)))/1000  
                             } else{
                               Dist <- NA
                             }
                             return(Dist)
                           })
      # Make a summary for mahalanobis distance estimates  
      MDminSumm <- c(MDmin,MDMinDistinKm)
      names(MDminSumm) <- c(names(MDmin),"DistinKm")
      writeRaster(MDminSumm,
                  paste0("./Results/Novelty/AllModels/",RCP,"/Seasonal/AllModels_",RCP,"_",YearUse,"_MDminSumm.tif"),
                  overwrite = TRUE)
      
      # Estimate for each Future Clime ensemble, where is the closest analogue using the Standarized Euclidean Distance
      SEDMin <- app(RCPFull,
                    function(TrgCellVals){#(TrgCellVals<-values(RCPFull)[1000,])
                      if(!is.na(TrgCellVals[1])){
                        SEDRast <- sum(((ClimNormMn-TrgCellVals)^2)/ClimNormSD)^0.5
                        out <- c(SED.Min = (minmax(SEDRast)['min',]), # Define the SEDmin value of a future cell to all Climate normal cell
                                 Cell = which(values(SEDRast)==(minmax(SEDRast)['min',]))) # Define  normal cell(s) that has(have) the SEDmin value
                        
                      } else{
                        out <- c(SED.Min = NA, # Define the SEDmin value of a future cell to all Climate normal cell
                                 Cell = NA) # Define  normal cell(s) that has(have) the SEDmin value
                      }
                      return(out)
                    })
      # Estimate for each Future Clime ensemble, the distance (in KM)  to the closest analogue using the Standarized Euclidean Distance
      CordsAll <- crds(RCPFull, df = TRUE,na.rm = FALSE) # get the coordinates
      SEDMinDistinKm <- rast(ClimNormMn[[1:3]]) # Make am empty SpatRaster File to summarize values
      values(SEDMinDistinKm) <- cbind(CordsAll,values(SEDMin[[2]])) # Add values (Coordinates) and the closest Cell
      names(SEDMinDistinKm) <- c("x","y","ID") # Renames so the distance function works (positions need to be x/y OR lat/lon)
      # Estimate the distance of each cell to the Closest analogue
      SEDMinDistinKm<- app(SEDMinDistinKm,
                           function(x){#(x<-values(SEDMinDistinKm)[1000,])
                             if(!is.na(x[3])){
                               Dist <- as.numeric(min(terra::distance(x = as.matrix(t(x[1:2])),
                                                                      y = as.matrix(CordsAll[x[3],]),
                                                                      lonlat = T)))/1000  
                             } else{
                               Dist <- NA
                             }
                             return(Dist)
                           })
      # Make a summary for SED distance estimates  
      SEDminSumm <- c(SEDMin,SEDMinDistinKm)
      names(SEDminSumm) <- c(names(SEDMin),"DistinKm")
      writeRaster(SEDminSumm,
                  paste0("./Results/Novelty/AllModels/",RCP,"/Seasonal/AllModels_",RCP,"_",YearUse,"_SEDminSumm.tif"),
                  overwrite = TRUE)      
    }
    rm(list = ls()[!ls()%in%c("Model","RCP","BioclimVars","BIOMES","ModelsAll","HistFiles1980_2005","MDtresh","SEDtresh",
                              "ClimNormMn","ClimNormSD","YearUse")]);gc()
  }
  gc()
}
