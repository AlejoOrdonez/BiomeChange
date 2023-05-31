rm(list=ls());gc()
require(terra)
require(snowfall)
require(maptools)
require(analogue)
setwd("/Users/au467796/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
# Map of the world
#data(wrld_simpl)
#wrld_simpl <- vect(wrld_simpl)
# Biomes Map
BIOMES <- vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")
BiomesNames <- data.frame(code = c(1:14,98,99,100),
                          Name = c("Tropical and Subtropical Moist Broadleaf Forests",
                                   "Tropical and Subtropical Dry Broadleaf Forests",
                                   "Tropical and Subtropical Coniferous Forests",
                                   "Temperate Broadleaf and Mixed Forests",
                                   "Temperate Coniferous Forests",
                                   "Boreal forests/Taiga",
                                   "Tropical and Subtropical Grasslands, Savannas, and Shrublands",
                                   "Temperate Grasslands, Savannas, and Shrublands",
                                   "Flooded Grasslands and Savannas",
                                   "Montane Grasslands and Shrublands",
                                   "Tundra",
                                   "Mediterranean Forests, Woodlands, and Scrub",
                                   "Deserts and Xeric Shrublands",
                                   "Mangroves",
                                   "Rock/ice",
                                   "Lake/River",
                                   "Ocean"))

setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
ModelsAll <- c("bcc-csm1-1", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0",
               "GISS-E2-H", "GISS-E2-R", "HadGEM2-ES", "MPI-ESM-LR", "NorESM1-M")


for (Model in c(ModelsAll)){#(Model <- ModelsAll[[1]])
  ClimNormList <- lapply(c("P.Seasonal","Tm.Seasonal"),
                         function(VarUse){#(VarUse <- "Tm.Seasonal")
                           # Load the 1960 to 2005 Historical data 
                           HistFiles1980_2005 <-lapply(paste0(Model,"_",1980:2005,"_",VarUse,".tif"),
                                                       function(x){#(x<- paste0(Model,"_",1960:2005,"_",VarUse,".tif")[1])
                                                         rast(paste0("./Data/CMIP5/Processed/Historical/Seasonal/",x))
                                                       })
                           # make the Climate normal file by adding the 2006 to 2010 Historical data (from each RCP) to HistFiles1980_2005
                           ClimNorm4RCP <- lapply(c("RCP26", "RCP46", "RCP60", "RCP85"),
                                                  function(RCP){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[2])
                                                    if(length(grep(Model,dir(paste0("./Data/CMIP5/Processed/",RCP,"/Seasonal"))))!=0){
                                                      # Load the 2006 to 2010 Historical data for an RCP
                                                      RCPList <- lapply(paste0(Model,"_",2006:2010,"_",VarUse,".tif"),
                                                                        function(x){#(x<- paste0(Model,"_",2006:2010,"_",VarUse,".tif")[1])
                                                                          if(x%in%dir(paste0("./Data/CMIP5/Processed/",RCP,"/Seasonal"))){
                                                                            out <- rast(paste0("./Data/CMIP5/Processed/",RCP,"/Seasonal/",x))
                                                                          } else{
                                                                            out <- NA
                                                                          }
                                                                          return(out)
                                                                        })
                                                      # Merge Historical and RCP lists
                                                      ClimNorm4RCPTmp <- c(HistFiles1980_2005,
                                                                           RCPList)
                                                      # Estimate the mean of each band for the Climate Normal
                                                      ClimNorm4RCPMn <- do.call("c",
                                                                                lapply(names(ClimNorm4RCPTmp[[1]]),
                                                                                       function(VarUse){#(VarUse <- names(ClimNorm4RCPTmp[[1]])[1])
                                                                                         app(do.call("c",
                                                                                                     lapply(ClimNorm4RCPTmp,
                                                                                                            function(x){x[[VarUse]]})),mean)
                                                                                       }))
                                                      # Estimate the SD of each band for the Climate Normal
                                                      ClimNorm4RCPSD <- do.call("c",
                                                                                lapply(names(ClimNorm4RCPTmp[[1]]),
                                                                                       function(VarUse){#(VarUse <- names(ClimNorm4RCPTmp[[1]])[1])
                                                                                         app(do.call("c",
                                                                                                     lapply(ClimNorm4RCPTmp,
                                                                                                            function(x){x[[VarUse]]})),sd)
                                                                                       }))
                                                      # Complie Mn and SD in one list
                                                      OutList <- list(RCP.Mn = ClimNorm4RCPMn,
                                                                      RCP.SD = ClimNorm4RCPSD)                           
                                                    } else {
                                                      OutList <- list(RCP.Mn = NA,
                                                                      RCP.SD = NA)
                                                    }
                                                    # Return the List value
                                                    return(OutList)
                                                  })
                           names(ClimNorm4RCP) <- c("RCP26", "RCP46", "RCP60", "RCP85")
                           return(ClimNorm4RCP)
                         })
  names(ClimNormList) <- c("P.Seasonal","Tm.Seasonal")
  
  # Estimate Novelty based on Future to - CLimateNormal distance
  Novelty <- lapply(c("RCP26", "RCP46", "RCP60", "RCP85"),
                    function(RCP){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[1])
                      ClimNormMn <- c(ClimNormList[[1]][[RCP]][["RCP.Mn"]],ClimNormList[[2]][[RCP]][["RCP.Mn"]])
                      ClimNormSD <- c(ClimNormList[[1]][[RCP]][["RCP.SD"]],ClimNormList[[2]][[RCP]][["RCP.SD"]])
                      BiomeBsLn <- extract(BIOMES[,"BIOME"],
                                           crds(ClimNormMn[[1]],df=T))
                      MDtreshDist <- apply(values(ClimNormMn),
                                           1,
                                           function(TrgCellVals){#(TrgCellVals<-values(ClimNormMn)[1,])
                                             (mahalanobis(values(ClimNormMn), TrgCellVals, cov(values(ClimNormMn))))
                                           })
                      MDtresh <- roc(object = MDtreshDist,
                                     groups = BiomeBsLn$BIOME)
                      rm(MDtreshDist);gc()
                      
                      SEDtreshDist <- apply(values(ClimNormMn),
                                            1,
                                            function(TrgCellVals){#(TrgCellVals<-values(ClimNormMn)[1,])
                                              as.numeric(values(sum(((ClimNormMn-TrgCellVals)^2)/ClimNormSD)^0.5))
                                            })
                      SEDtresh <- roc(object = SEDtreshDist,
                                      groups = BiomeBsLn$BIOME)
                      rm(SEDtreshDist);gc()
                      gc()
                      # Estmate the distance to each period
                      NoveltyPerPer <- lapply(seq(2099,2299,by=50),
                                              function(YearUse){#(YearUse <- seq(2099,2299,by=50)[5])
                                                # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
                                                RCPFullList <- lapply(c("P.Seasonal","Tm.Seasonal"),
                                                                      function(VarUse){#(VarUse <- "Tm.Seasonal")
                                                                        RCPList <- lapply(paste0(Model,"_",c(c(YearUse-9):YearUse),"_",VarUse,".tif"),
                                                                                          function(x){#(x<- paste0(Model,"_",2006:2010,"_",VarUse,".tif")[1])
                                                                                            if(x%in%dir(paste0("./Data/CMIP5/Processed/",RCP,"/Seasonal"))){
                                                                                              out <- rast(paste0("./Data/CMIP5/Processed/",RCP,"/Seasonal/",x))
                                                                                            } else{
                                                                                              out <- NA
                                                                                            }
                                                                                            return(out)
                                                                                          })
                                                                        # Estimate the mean of each band for the 
                                                                        RCPMnTmp <- Reduce("+",RCPList)/length(RCPList)
                                                                      })
                                                RCPFullList
                                                # Turn the RCP List into a SpatRaster 
                                                RCPFull <- do.call("c",RCPFullList)
                                                # Estimate for each Future Clime ensemble, where is the closest analogue using the mahalanobis Distance
                                                MDmin <- app(RCPFull,
                                                             function(TrgCellVals){
                                                               MDAll <- (mahalanobis(values(ClimNormMn), TrgCellVals, cov(values(ClimNormMn))))
                                                               c(MD.min=min(MDAll), # Define the MDmin value of a future cell to all Climate normal cell
                                                                 Cell=which(MDAll==min(MDAll)))
                                                             })
                                                # Estimate for each Future Clime ensemble, the distance (in KM)  to the closest analogue using the mahalanobis Distance
                                                CordsAll <- crds(RCPFull, df=T) # get the coordinates
                                                MDMinDistinKm <- rast(ClimNormMn[[1:3]]) # Make am empty SpatRaster File to summarize values
                                                values(MDMinDistinKm) <- cbind(CordsAll,values(MDmin[[2]])) # Add values (Coordinates) and the closest Cell
                                                names(MDMinDistinKm) <- c("x","y","ID") # Renames so the distance function works (positions need to be x/y OR lat/lon)
                                                # Estimate the distance of each cell to the Closest analogue
                                                MDMinDistinKm<- app(MDMinDistinKm,
                                                                    function(x){
                                                                      as.numeric(min(distance(x = as.matrix(t(x[1:2])),
                                                                                              y = as.matrix(CordsAll[x[3],]),
                                                                                              lonlat = T)))/1000
                                                                    })
                                                # Make a summary for mahalanobis distance estimates  
                                                MDminSumm <- c(MDmin,MDMinDistinKm)
                                                names(MDminSumm) <- c(names(MDmin),"DistinKm")
                                                
                                                # Estimate for each Future Clime ensemble, where is the closest analogue using the Standarized Euclidean Distance
                                                SEDMin <- app(RCPFull,
                                                              function(TrgCellVals){
                                                                SEDRast <- sum(((ClimNormMn-TrgCellVals)^2)/ClimNormSD)^0.5
                                                                c(SED.Min = (minmax(SEDRast)['min',]), # Define the SEDmin value of a future cell to all Climate normal cell
                                                                  Cell = which(values(SEDRast)==(minmax(SEDRast)['min',])) # Define  normal cell(s) that has(have) the SEDmin value
                                                                )
                                                              })
                                                # Estimate for each Future Clime ensemble, the distance (in KM)  to the closest analogue using the Standarized Euclidean Distance
                                                CordsAll <- crds(RCPFull, df=T)# get the coordinates
                                                SEDMinDistinKm <- rast(ClimNormMn[[1:3]])  # Make am empty SpatRaster File to summarize values
                                                values(SEDMinDistinKm) <- cbind(CordsAll,values(SEDMin[[2]])) # Add values (Coordinates) and the closest Cell
                                                names(SEDMinDistinKm) <- c("x","y","ID") # Renames so the distance function works (positions need to be x/y OR lat/lon)
                                                # Estimate the distance of each cell to the Closest analogue
                                                SEDMinDistinKm<- app(SEDMinDistinKm,
                                                                     function(x){
                                                                       as.numeric(min(distance(x = as.matrix(t(x[1:2])),
                                                                                               y = as.matrix(CordsAll[x[3],]),
                                                                                               lonlat = T)))/1000
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
                                                                                SEDtresh$roc$Combined$optimal))
                                                #Save the Output
                                                saveRDS(Out.List,
                                                        paste0("./Results/Novelty/Seasonal/",Model,"_",RCP,"_",YearUse,".rds"))
                                                #~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Results/Novelty/bcc-csm1-1_RCP26_2099.rds
                                                return(Novelty)
                                              })
                    })
  rm(ClimNormList);gc()
}