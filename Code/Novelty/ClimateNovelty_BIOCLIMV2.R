rm(list=ls());gc()
require(terra)
require(snowfall)
require(maptools)
require(analogue)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")

# BioclimVars to use
BioclimVars <- c(5, # Max Temperature of Warmest Month
                 6, # Min Temperature of Coldest Month
                 13,# Precipitation of Wettest Month
                 14)
# ,# Precipitation of Driest Month
#                  8, # Mean Temperature of Wettest Quarter
#                  9, # Mean Temperature of Driest Quarter
#                  10,# Mean Temperature of Warmest Quarter
#                  11,# Mean Temperature of Coldest Quarter
#                  16,# Precipitation of Wettest Quarter
#                  17,# Precipitation of Driest Quarter
#                  18,# Precipitation of Warmest Quarter
#                  19)# Precipitation of Coldest Quarter

# Biomes Map
BIOMES <- vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")

ModelsAll <- c("bcc-csm1-1", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0",
               "GISS-E2-H", "GISS-E2-R", "HadGEM2-ES", "MPI-ESM-LR", "NorESM1-M")


for (Model in c(ModelsAll)){#(Model <- ModelsAll[[1]])
  # Load the 1960 to 2005 Historical data 
  HistFiles1980_2005 <-lapply(paste0(Model,"_",1980:2005,"_BIOCLIM.tif"),
                              function(x){#(x<- paste0(Model,"_",1960:2005,"_BIOCLIM.tif")[1])
                                rast(paste0("./Data/CMIP5/Processed/Historical/BIOCLIM/",x))[[BioclimVars]]
                              })
  
  # Estimate Novelty based on Future to - CLimate Normal distance
  for (RCP in c("RCP26", "RCP46", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[1])
    if(length(grep(Model,dir(paste0("./Data/CMIP5/Processed/",RCP,"/BIOCLIM"))))!=0){
      # Load the 2006 to 2010 Historical data for an RCP to  Create a climate normal raster for each evaluated variable  
      RCPList <- lapply(paste0(Model,"_",2006:2010,"_BIOCLIM.tif"),
                        function(x){#(x<- paste0(Model,"_",2006:2010,"_BIOCLIM.tif")[1])
                          rast(paste0("./Data/CMIP5/Processed/",RCP,"/BIOCLIM/",x))[[BioclimVars]]
                        })
      # Merge Historical and RCP lists
      ClimNorm4RCPTmp <- c(HistFiles1980_2005,
                           RCPList)
      # Estimate the mean of each band for the Climate Normal
      ClimNormMn <- do.call("c",
                            lapply(names(ClimNorm4RCPTmp[[1]]),
                                   function(VarUse){#(VarUse <- names(ClimNorm4RCPTmp[[1]])[1])
                                     app(do.call("c",
                                                 lapply(ClimNorm4RCPTmp,
                                                        function(x){x[[VarUse]]})),mean)
                                   }))
      names(ClimNormMn) <- paste0("bio",BioclimVars)
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
      names(ClimNormSD) <- paste0("bio",BioclimVars)
      #Crop oceans
      ClimNormSD <- mask(ClimNormSD,BIOMES[BIOMES$BIOME<15,])
      
      # Estimate Novelty threshold based on Climate-Normal distance                                                   
      # Define the "BIOME" for each location
      BiomeBsLn <- extract(BIOMES[,"BIOME"],
                           crds(ClimNormMn[[1]],df=T))
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
      # Estimate pairwise SED differences for the Climate-Normal period
      SEDtreshDist <- apply(values(ClimNormMn),
                            1,
                            function(TrgCellVals){#(TrgCellVals<-values(ClimNormMn)[1,])
                              as.numeric(values(sum(((ClimNormMn-TrgCellVals)^2)/ClimNormSD)^0.5))
                            })
      # Estimate the SED distance based no-anlaogue threshold 
      SEDtresh <- roc(object = SEDtreshDist,
                      groups = BiomeBsLn$BIOME)
      rm(SEDtreshDist);gc()
      
      
      # Estimate Novelty by comparing future climate to the Climate-Normal
      # Estimate the distance to each period
      rm(list = ls()[!ls()%in%c("Model","RCP","BioclimVars","BIOMES","ModelsAll","HistFiles1980_2005","MDtresh","SEDtresh",
                                "ClimNormMn","ClimNormSD")]);gc()
      for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[2])
        # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
        RCPList <- lapply(paste0(Model,"_",c(c(YearUse-9):YearUse),"_BIOCLIM.tif"),
                          function(x){#(x<- paste0(Model,"_",2006:2010,"_",VarUse,".tif")[1])
                            rast(paste0("./Data/CMIP5/Processed/",RCP,"/BIOCLIM/",x))[[BioclimVars]]})
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
        
        # final summary all Values
        Out.List <- list(RCP = RCP,
                         Year = YearUse,
                         MDSumm = list(MDminSumm,
                                       Tresh = MDtresh$roc$Combined$optimal),
                         SEDSumm = list(SEDminSumm,
                                        Tresh = SEDtresh$roc$Combined$optimal))
        #Save the Output
        saveRDS(Out.List,
                paste0("./Results/Novelty/BIOCLIM/",Model,"_",RCP,"_",YearUse,".rds"))
        rm(list = ls()[!ls()%in%c("Model","RCP","BioclimVars","BIOMES","ModelsAll","HistFiles1980_2005","MDtresh","SEDtresh",
                                  "ClimNormMn","ClimNormSD","YearUse")]);gc()
      }
    }
  }
}



