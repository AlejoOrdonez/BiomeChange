rm(list=ls());gc()
require(terra)
require(snowfall)
require(maptools)
data(wrld_simpl)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
ModelsAll <- c("bcc-csm1-1", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0",
               "GISS-E2-H", "GISS-E2-R", "HadGEM2-ES", "MPI-ESM-LR", "NorESM1-M")
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/Data/CimateData/CMIP5")

# for (Model in c(ModelsAll)){#(Model <- ModelsAll[[1]])
(Model <- ModelsAll[[1]])
  CLimNormList <- lapply(c("P.Seasonal","Tm.Seasonal"),
                         function(VarUse){#(VarUse <- "Tm.Seasonal")
                           # Load the 1960 to 2005 Historical data 
                           HistFiles1980_2005 <-lapply(paste0(Model,"_",1980:2005,"_",VarUse,".tif"),
                                                       function(x){#(x<- paste0(Model,"_",1960:2005,"_",VarUse,".tif")[1])
                                                         rast(paste0("./Processed/Historical/Seasonal/",x))
                                                       })
                           # make the CLimate normal file by adding the 2006 to 2010 Historical data (from each RCP) to HistFiles1980_2005
                           ClimNorm4RCP <- lapply(c("RCP26", "RCP46", "RCP60", "RCP85"),
                                                  function(RCP){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[2])
                                                    if(length(grep(Model,dir(paste0("./Processed/",RCP,"/Seasonal"))))!=0){
                                                      # Load the 2006 to 2010 Historical data for an RCP
                                                      RCPList <- lapply(paste0(Model,"_",2006:2010,"_",VarUse,".tif"),
                                                                        function(x){#(x<- paste0(Model,"_",2006:2010,"_",VarUse,".tif")[1])
                                                                          if(x%in%dir(paste0("./Processed/",RCP,"/Seasonal"))){
                                                                            out <- rast(paste0("./Processed/",RCP,"/Seasonal/",x))
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

  



# Novelty <- lapply(c("RCP26", "RCP46", "RCP60", "RCP85"),
#                   function(RCP){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[1])
  (RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[1])
                      ClimNormMn <- c(CLimNormList[[1]][[RCP]][["RCP.Mn"]],CLimNormList[[2]][[RCP]][["RCP.Mn"]])
                      ClimNormSD <- c(CLimNormList[[1]][[RCP]][["RCP.SD"]],CLimNormList[[2]][[RCP]][["RCP.SD"]])
#                     lapply(seq(2099,2299,by=50),
#                            function(YearUse){#(YearUse <- seq(2099,2299,by=50)[5])
#                            })
#                   })
  (YearUse <- seq(2099,2299,by=50)[5])
                             a <- Sys.time()

                             # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
                             RCPFullList <- lapply(c("P.Seasonal","Tm.Seasonal"),
                                                   function(VarUse){#(VarUse <- "Tm.Seasonal")
                                                     RCPList <- lapply(paste0(Model,"_",c(c(YearUse-9):YearUse),"_",VarUse,".tif"),
                                                                       function(x){#(x<- paste0(Model,"_",2006:2010,"_",VarUse,".tif")[1])
                                                                         if(x%in%dir(paste0("./Processed/",RCP,"/Seasonal"))){
                                                                           out <- rast(paste0("./Processed/",RCP,"/Seasonal/",x))
                                                                         } else{
                                                                           out <- NA
                                                                         }
                                                                         return(out)
                                                                       })
                                                     # Estimate the mean of each band for the 
                                                     RCPMnTmp <- Reduce("+",RCPList)/length(RCPList)
                                                   })
                             RCPFull <- do.call("c",RCPFullList)
                             # Estimate for each Future clime ensemble, where is the closest analogue using the mahalanobis Distance
                             b <- Sys.time()
                             MDmin <- app(RCPFull,
                                       function(TrgCellVals){
                                         MDAll <- (mahalanobis(values(ClimNormMn), TrgCellVals, cov(values(ClimNormMn))))
                                         c(MD.min=min(MDAll), # Define the MDmin value of a future cell to all Climate normal cell
                                           Cell=which(MDAll==min(MDAll)))
                                       })
                             Sys.time() - b
                             plot(MDmin[[1]])
                             plot(wrld_simpl,add=T)
                             hist(MDmin[[1]])
                             # Estimate for each Future clime ensemble, the distance (in KM)  to the closest analogue using the mahalanobis Distance
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
                             plot(MDMinDistinKm)
                             plot(wrld_simpl,add=T)
                             
                             
                             # Estimate for each Future clime ensemble, where is the closest analogue using the Standarized Euclidean Distance
                             RCPMnTmpDatFrm <- values(RCPFull,dataframe=T) # Make the Future Data a DataFrame
                             b <- Sys.time()
                             SEDMin <- app(RCPFull,
                                           function(TrgCellVals){
                                             SEDRast <- sum(((ClimNormMn-TrgCellVals)^2)/ClimNormSD)^0.5
                                             c(SED.Min = (minmax(SEDRast)['min',]), # Define the SEDmin value of a future cell to all Climate normal cell
                                               Cell = which(values(SEDRast)==(minmax(SEDRast)['min',])) # Define  normal cell(s) that has(have) the SEDmin value
                                               )
                                           })
                             Sys.time() - b
                             plot(SEDMin[[1]])
                             plot(wrld_simpl,add=T)
                             hist(SEDMin[[1]])
                             
                             # Estimate for each Future clime ensemble, the distance (in KM)  to the closest analogue using the Standarized Euclidean Distance
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
                             plot(SEDMinDistinKm)
                             plot(wrld_simpl,add=T)
                             
