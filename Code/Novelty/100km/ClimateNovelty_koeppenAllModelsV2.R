rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")

# Biomes Map
BIOMES <- vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")
BIOMES <- project(BIOMES,"+proj=eck4")

# Create a 100 x100 km raster
WGSRast <- rast(nrows=180, ncols=360) 
eck4Rast <- project(WGSRast,"+proj=eck4")
eck4Rast <- rast(extent=ext(eck4Rast),resolution=100000,crs="+proj=eck4")
rm(WGSRast);gc()

# Estimate Novelty based on Future to - CLimate Normal distance
#####
#####
  # Estimate the baseline mean and SD of climate variables
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[2])
  # Load the 1980 to 2010 Historical data for an RCP to Create a climate normal raster for each evaluated variable  
  ClimNorm4RCPTmp <- lapply(1980:2010,
                            function(YearUse){#(YearUse <- c(1980:2010)[1])
                            tmp <- rast(paste0("./Data/CMIP5/Processed/AllModelSumm1Arcmin/",RCP,"/koeppen_geiger/AllModels_",RCP,"_",YearUse,".tif"))
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
  BiomeBsLn <- rasterize(BIOMES,ClimNormMn[[1]],field = "BIOME")
  BiomeBsLn <- app(BiomeBsLn,function(x){ifelse(is.na(x),NA,ifelse(x>14,NA,x))})
  #Crop oceans
  ClimNormMn <- mask(ClimNormMn,BiomeBsLn)
  ClimNormSD <- mask(ClimNormSD,BiomeBsLn)
#####
#####
# Estimate Novelty threshold based on Climate-Normal distance
  if(!paste0("AllModels_",RCP,"_TreshSumm.rds")%in%dir(paste0("./Results2/Novelty/AllModels/",RCP,"/koeppen_geiger/"))){
    TimeIn <- Sys.time()
    # Estimate the mahalanobis distance based no-anlaogue threshold - use a resample approch where all Biomes have the same number of Observations
    MDtresh <- lapply(1:1000,
                      function(x){
                        BiomeSamp <- spatSample(BiomeBsLn, 40, "stratified", as.points=F,cells=FALSE)
                        MDtreshDist <- apply((ClimNormMn[][BiomeSamp$cell,]),
                                             1,
                                             function(TrgCellVals){#(TrgCellVals<-ClimMn[1,])
                                               out <- mahalanobis((ClimNormMn[][BiomeSamp$cell,]),
                                                                  TrgCellVals,
                                                                  cov((ClimNormMn[][BiomeSamp$cell,])))
                                               out <- round(out,3)
                                               return(out)
                                             })
                        MDtresh <- roc(object = MDtreshDist,
                                       groups = BiomeSamp$lyr.1)
                        MDtresh$roc$Combined$optimal
                        
                      })
    rm(MDtreshDist);gc()
     1.379Sys.time()-TimeIn
    # Estimate pairwise SED differences for the Climate-Normal period
    TimeIn <- Sys.time()
    SEDtreshDist <- apply(ClimMn,
                          1,
                          function(TrgCellVals){#(TrgCellVals<-ClimMn[1,])
                            out <- as.numeric(values(sum(((ClimNormMn-TrgCellVals)^2)/ClimNormSD)^0.5,na.rm=T))
                            out <- round(out,3)
                            return(out)
                          })
    # Estimate the SED distance based no-anlaogue threshold 
    
    
    SEDtresh <- roc(object = SEDtreshDist,
                    groups = values(BiomeBsLn, na.rm=T))
    rm(SEDtreshDist);gc()
    Sys.time()-TimeIn
    # final summary all Values
    Out.List <- list(RCP = RCP,
                     MDSummTresh = MDtresh$roc$Combined$optimal,
                     SEDSummTresh = SEDtresh$roc$Combined$optimal)
    #Save the Output
    saveRDS(Out.List,
            paste0("./Results2/Novelty/AllModels/",RCP,"/koeppen_geiger/AllModels_",RCP,"_TreshSumm.rds"))
  }
  rm(list = ls()[!ls()%in%c("RCP","BIOMES",
                            "MDtresh","SEDtresh","Out.List",
                            "ClimNormMn","ClimNormSD")]);gc()
#####
#####
  # Estimate Novelty by comparing future climate to the Climate-Normal
  for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[2])
    if(!paste0("AllModels_",RCP,"_",YearUse,"_SEDminSumm.tif")%in%dir(paste0("./Results2/Novelty/AllModels/",RCP,"/koeppen_geiger/"))){
      # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
      RCPList <- lapply(c(YearUse-9):YearUse,
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
                  paste0("./Results2/Novelty/AllModels/",RCP,"/koeppen_geiger/AllModels_",RCP,"_",YearUse,"_MDminSumm.tif"),
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
                  paste0("./Results2/Novelty/AllModels/",RCP,"/koeppen_geiger/AllModels_",RCP,"_",YearUse,"_SEDminSumm.tif"),
                  overwrite = TRUE)      
    }
    rm(list = ls()[!ls()%in%c("Model","RCP","BioclimVars","BIOMES","ModelsAll","HistFiles1980_2005","MDtresh","SEDtresh",
                              "ClimNormMn","ClimNormSD","YearUse")]);gc()
  }
  gc()
}

