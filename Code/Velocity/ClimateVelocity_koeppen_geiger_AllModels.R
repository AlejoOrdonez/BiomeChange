rm(list=ls());gc()
require(terra)
require(snowfall)
require(maptools)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
source("./Code/Velocity/VelocityFnc.R")

# Biomes Map
BIOMES <- vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")



# Estimate Novelty based on Future to - CLimate Normal distance
for (RCP in c("RCP26", "RCP46", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[1])
  # Load the 1980 to 2010 Historical data for an RCP to  Create a climate normal raster for each evaluated variable  
  RCPList <- lapply(c(1980:2009),
                    function(YearUse){#(YearUse <- c(1980:2009)[1])
                      rast(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/koeppen_geiger/AllModels_",RCP,"_",YearUse,".tif"))
                    })
  # Estimate the mean of each band for the Climate Normal
  ClimNormMn <- do.call("c",
                        lapply(names(RCPList[[1]]),
                               function(VarUse){#(VarUse <- names(RCPList[[1]])[1])
                                 x <- do.call("c",
                                              lapply(RCPList,
                                                     function(x){x[[VarUse]]}))
                                 app(x,
                                     fun=function(i, ff) ff(i),
                                     cores = 3, ff= function(i) {mean(i)})
                               }))
  names(ClimNormMn) <- names(RCPList[[1]])
  #Crop oceans
  ClimNormMn <- mask(ClimNormMn,BIOMES[BIOMES$BIOME<15,])
  #Clean the Memory
  rm(RCPList);gc()

  for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[3])
    ### Velocity per variable
    for(var in names(ClimNormMn)){#(var <- names(ClimNormMn)[1])
      # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
      RCPList <- lapply(c(1980:YearUse),
                        function(x){#(x <- c(1980:YearUse)[1])
                          rast(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/koeppen_geiger/AllModels_",RCP,"_",x,".tif"),
                               lyrs = var )
                        })
      a <- crds(RCPList[[1]])
      d <- distance(a,lonlat=TRUE)
      # Step 0: create a raster with the time series for a given variable 
      TimeSerRast <- do.call("c",RCPList)
      TimeSerRast2 <- mask(TimeSerRast,BIOMES[BIOMES$BIOME<15,])
      # Step 1a: Estimate the Temporal gradient as the slope of the time series
      TimeHetARM <- app(TimeSerRast,
                        fun=function(i, ff) ff(i),
                        cores = 3,
                        ff = function(x){
                          if(!is.na(x[1]) & length(unique(x))!=1){
                            tmpDtaFrm <- data.frame(prop = x,
                                                    Time = c(1:length(x)))
                            TimMod <- nlme::gls(prop~Time, data = tmpDtaFrm,
                                                correlation = nlme::corARMA(p=1), 
                                                method ="ML",
                                                control=list(opt = "optim",
                                                             msMaxIter = 100000))
                            coef(summary(TimMod))["Time",c("Value")]
                          }
                          if(length(unique(x))==1){
                          0
                            }
                          else{
                            NA
                          }
                        })
      TimeHetARM <- mask(TimeHetARM,BIOMES[BIOMES$BIOME<15,])
      #  Step 1b: Estimate the Temporal gradient as the Median of inter-anual changes        
      TimeHetIntAnn<- app(TimeSerRast,
                          fun=function(i, ff) ff(i,method = "Anomaly1"),
                          cores = 3,
                          ff = TempGradFnc)
      TimeHetIntAnn <-  mask(TimeHetIntAnn,BIOMES[BIOMES$BIOME<15,])
      # Step 1c: Estimate the Temporal gradient as the anomaly divided by time
      TimeHetAnom <- abs(TimeSerRast[[dim(TimeSerRast2)[3]]]-TimeSerRast[[1]])/dim(TimeSerRast2)[3]
      TimeHetAnom <-  mask(TimeHetAnom,BIOMES[BIOMES$BIOME<15,])
      
      # Step 2.: Estimate the Spatial gradient -  Here we use the ClimNormMn raster
      SpatHetTmp <- SpatHetFnc(ClimNormMn[[var]],Dist=28) # NEED TO DEFINE THE DISTACNE BETWEEN CELLS 
      
      # Step 3: Estimate the Velocity as the ratio between spatial and temporal gradients 
      # Step 3a: estimate the Velocity using a slope of the time series derived Temporal gradient
      VelocityARM <- VelocityFnc(TimeHetARM,
                                 SpatHetTmp)
      # Step 3b: estimate the Velocity using a sMedian of inter-anual changes derived Temporal gradient
      VelocityIntAnn <- VelocityFnc(TimeHetIntAnn,
                                    SpatHetTmp)
      # Step 3c: estimate the Velocity using the Start-End anomaly derived Temporal gradient
      VelocityAnom <- VelocityFnc(TimeHetAnom,
                                  OutRast[['SpatHet']])
      
      # Step 4.: Estimate the Bearing using the degrees from north of the vector-sum used to estimate the Spatial gradient
      BearingTmp <- BearingFnc(c(ClimNormMn[[var]],
                                 RCPList[[length(RCPList)]][[var]]))
      # Step 5.: Create a summary for a variable
      OutRast <- c(TimeHetARM,TimeHetIntAnn,TimeHetAnom,
                   SpatHetTmp,BearingTmp,
                   VelocityARM,VelocityIntAnn,
                   VelocityAnom)
      names(OutRast) <- c("TimeH.ARM", "TimeH.IntAnn", "TimeH.StEnAnn", "SpatHet", "Bearing", "Vel.ARM", "Vel.IntAnn", "Vel.StEnAnn")
      #Save the Output
      writeRaster(OutRast,
                  filename = paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_1980to",YearUse,".tif"),
                  overwrite = TRUE)
      
      }
  }
        
        ### Estimate multivariate metrics - Displacement
        # Step 1a. Load the ARM derived velocities        
        Vel.ARM.List <- lapply(names(ClimNormMn),
                               function(var){#(var <- names(ClimNormMn)[1])
                                 rast(paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_1980to",YearUse,".tif"),
                                      lyrs = "Vel.ARM")
                               })
        # Step 2a. Estimate the displacement as the median velocity  
        Displ.ARM <- 10^app(log10(do.call("c",Vel.ARM.List)),
                            fun=function(i, ff) ff(i),
                            cores =3,
                            ff=function(x){median(x[x > -2 & x < 2])})
        
        # Step 1b. Load the Inter annual change derived velocities        
        Vel.IntAnn.List <- lapply(names(ClimNormMn),
                                  function(var){#(var <- names(ClimNormMn)[1])
                                    rast(paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_1980to",YearUse,".tif"),
                                         lyrs = "Vel.IntAnn")
                                  })
        # Step 2b. Estimate the displacement as the median velocity  
        Displ.IntAnn <- 10^app(log10(do.call("c",Vel.IntAnn.List)),
                               fun=function(i, ff) ff(i),
                               cores =3,
                               ff=function(x){median(x[x > -2 & x < 2])})
        
        # Step 1c. Load the Start-End anomaly derived velocities        
        Vel.StEnAnn.List <- lapply(names(ClimNormMn),
                                   function(var){#(var <- names(ClimNormMn)[1])
                                     abs(rast(paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_1980to",YearUse,".tif"),
                                              lyrs = "Vel.StEnAnn"))
                                   })
        # Step 2c. Estimate the displacement as the median velocity  
        Displ.StEnAnn <- 10^app(log10(do.call("c",Vel.StEnAnn.List)),
                                fun=function(i, ff) ff(i),
                                cores =3,
                                ff=function(x){median(x[x > -2 & x < 2])})
        
        ## Estimate multivariate metrics- Divergence
        # Step 1.: Load the bearings
        BearingList <- lapply(names(ClimNormMn),
                              function(var){#(var <- names(ClimNormMn)[1])
                                rast(paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_1980to",YearUse,".tif"),
                                     lyrs = "Bearing")
                              })
        # Estimate the Divergence as the median difference between bearings -  median angle between bearings
        Divergence <- app(do.call("c",BearingList),
                          function(x){
                            if(is.na(x[1])){out<-NA
                            } else{x <- x[x!=361]
                            #x <- round(x,1) 
                            y <- dist(na.omit(x))
                            y[y>180] <- 360 - y[y>180]
                            out <- median(y) }
                            return((out))
                          })
        # Save the Multivariate Output
        OutRastMultVar <- c(Displ.ARM, Displ.IntAnn,Displ.StEnAnn,Divergence)
        names(OutRastMultVar) <- c("Displ.ARM","Displ.IntAnn","Displ.StEnAnn","Dive")
        
        ## Save the Output
        writeRaster(OutRastMultVar,
                    filename = paste0("./Results/Displacement_Divergence/koeppen_geiger/",Model,"_DispDiv_",RCP,"_1980to",YearUse,".tif"),
                    overwrite = TRUE)
        ## Clean the Memory
        rm(list = ls()[!ls()%in%c("Model","RCP","koeppen_geigerVars","BIOMES","ModelsAll",
                                  "HistFiles1980_2005","ClimNormMn",
                                  "TempGradFnc", "SpatHetFnc","VelocityFnc","BearingFnc")]);gc()
      }
    }
  }
}



