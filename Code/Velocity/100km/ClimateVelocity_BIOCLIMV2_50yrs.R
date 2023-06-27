rm(list=ls());gc()
require(terra)
require(snowfall)
require(maptools)
setwd("/Users/alejandroordonez/Downloads/BiomeChange/")
source("./Code/Velocity/VelocityFnc.R")
# BioclimVars to use
BioclimVars <- c(5, # Max Temperature of Warmest Month
                 6, # Min Temperature of Coldest Month
                 13,# Precipitation of Wettest Month
                 14) # Precipitation of Driest Month

# Biomes Map
BIOMES <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_100km.tif")

# Estimate Novelty based on Future to - Climate Normal distance
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[2])
  for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[1])
    ### Velocity per variable
    for(var in paste0("bio",BioclimVars)){#(var <- paste0("bio",BioclimVars)[1])
      if(!paste0("AllModels_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif")%in% dir(paste0("./Results2/Velocity/AllModels/",RCP,"/BIOCLIM/"))){
        # Step 0: create a raster with the time series for a given variable 
        # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
        RCPList <- lapply(paste0("AllModels_",RCP,"_",c(YearUse-50):c(YearUse),".tif"),
                          function(x){#(x<- paste0("AllModels_",RCP,"_",c(YearUse-50):c(YearUse),".tif")[2])
                            tmp <- rast(paste0("./Data/CMIP5/Processed/",RCP,"/BIOCLIM/",x),
                                        lyrs = as.numeric(gsub("bio","",var)))
                            # Project and match the target resolution
                            tmp <- project(tmp,"+proj=eck4")
                            tmp<- resample(tmp,
                                           BIOMES,
                                           method = "near")
                            #Crop oceans
                            tmp <- mask(tmp,BIOMES)
                            
                            return(tmp)
                          })
        TimeSerRast <- do.call("c",RCPList)
        names(TimeSerRast) <- c(YearUse-50):c(YearUse)
        rm(RCPList);gc()
        # Step 1a: Estimate the Temporal gradient as the slope of the time series
        TimeHetARM <- app(TimeSerRast,
                          fun=function(i, ff) ff(i),
                          cores = 3,
                          ff = function(x){
                            if(!is.na(x[1]) & length(unique(x))>1){
                              tmpDtaFrm <- data.frame(prop = x,
                                                      Time = c(1:length(x)))
                              TimMod <- nlme::gls(prop~Time, data = tmpDtaFrm,
                                                  correlation = nlme::corARMA(p=1), 
                                                  method ="ML",
                                                  control=list(opt = "optim",
                                                               msMaxIter = 100000))
                              coef(summary(TimMod))["Time",c("Value")]
                            } else{
                              ifelse(length(unique(x))==1 & !is.na(x[1]),
                                     0,
                                     NA)}
                          })
        
        #  Step 1d: Estimate the Temporal gradient as the Median of inter-annual changes        
        TimeHetIntAnn<- app(TimeSerRast,
                            fun=function(i, ff) ff(i),
                            cores = 3,
                            ff = function(x){
                              if(!is.na(x[1]) & length(unique(x))>1){
                                TimMod <- x[-1] - x[-length(x)]
                                Out <- median(TimMod)
                              } else{
                                ifelse(length(unique(x))==1 & !is.na(x[1]),
                                       0,
                                       NA)}
                            })
        #  Step 1c: Estimate the Temporal gradient as the Dif between Start and end conditions        
        #Create a climate normal raster for each evaluated variable  based on the data for the 30 years before the last start period 
        ClimNormMn <- lapply(paste0("AllModels_",RCP,"_",c(YearUse-80):c(YearUse-50),".tif"),
                             function(x){#(x<- paste0("AllModels_",RCP,"_",c(YearUse-50):c(YearUse),".tif")[1])
                               tmp <- rast(paste0("./Data/CMIP5/Processed/",RCP,"/BIOCLIM/",x),
                                           lyrs = as.numeric(gsub("bio","",var)))
                               # Project and match the target resolution
                               tmp <- project(tmp,"+proj=eck4")
                               tmp<- resample(tmp,
                                              BIOMES,
                                              method = "near")
                               #Crop oceans
                               tmp <- mask(tmp,BIOMES)
                               return(tmp)
                             })
        ClimNormMn <- mean(do.call("c",ClimNormMn))
        
        #Create a Future climate raster for each evaluated variable  based on a 30 years period centered at the end of the period  
        ClimFut <- lapply(paste0("AllModels_",RCP,"_",c(YearUse-30):c(YearUse),".tif"),
                          function(x){#(x<- paste0("AllModels_",RCP,"_",c(YearUse-50):c(YearUse),".tif")[1])
                            tmp <- rast(paste0("./Data/CMIP5/Processed/",RCP,"/BIOCLIM/",x),
                                        lyrs = as.numeric(gsub("bio","",var)))
                            # Project and match the target resolution
                            tmp <- project(tmp,"+proj=eck4")
                            tmp<- resample(tmp,
                                           BIOMES,
                                           method = "near")
                            #Crop oceans
                            tmp <- mask(tmp,BIOMES)
                            return(tmp)
                          })
        ClimFut <- do.call("c",ClimFut)
        ClimFut <- mean(ClimFut)
        TimeHetAnn <- (ClimFut - ClimNormMn)/50
        # Step 2.: Estimate the Spatial gradient -  Here we use the ClimNormMn raster
        SpatHetTmp <- SpatHetFnc(RastIn = ClimNormMn,
                                 Dist = res(BIOMES)[1]/1000)
        # Step 3: Estimate the Velocity as the ratio between spatial and temporal gradients 
        # Step 3a: estimate the Velocity using a slope of the time series derived Temporal gradient
        VelocityARM <- VelocityFnc(TimeHetARM,
                                   SpatHetTmp)
        # Step 3b: estimate the Velocity using a Median of inter-anual changes derived Temporal gradient
        VelocityIntAnn <- VelocityFnc(TimeHetIntAnn,
                                      SpatHetTmp)
        # Step 3cc: estimate the Velocity using a sStrat end Temporal gradient
        VelocityAnn <- VelocityFnc(TimeHetAnn,
                                   SpatHetTmp)
        # Step 4.: Estimate the Bearing using the degrees from north of the vector-sum used to estimate the Spatial gradient
        BearingTmp <- BearingFnc(c(ClimNormMn,
                                   ClimFut))
        # Step 5.: Create a summary for a variable
        OutRast <- c(TimeHetARM, TimeHetIntAnn, TimeHetAnn,
                     SpatHetTmp, BearingTmp,
                     VelocityARM, VelocityIntAnn, VelocityAnn)
        names(OutRast) <- c("TimeH.ARM", "TimeH.IntAnn","TimeH.Ann",
                            "SpatHet", "Bearing",
                            "Vel.ARM", "Vel.IntAnn", "Vel.Ann")
        #Save the Output
        writeRaster(OutRast,
                    filename = paste0("./Results2/Velocity/AllModels/",RCP,"/BIOCLIM/","AllModels_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
                    overwrite = TRUE)
        
        # clean memory
        rm(list = c("TimeSerRast",
                    "TimeHetARM","TimeHetIntAnn","TimeHetAnn",
                    "ClimNormMn","ClimFut",
                    "SpatHetTmp",
                    "VelocityARM","VelocityIntAnn","VelocityAnn",
                    "BearingTmp",
                    "OutRast"))
        gc()                    
      }
    }
  }
}



for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[4])
  for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[1])
    ### Estimate multivariate metrics - Displacement
    if(!paste0("AllModels_DispDiv_",RCP,"_",c(YearUse-50),"to",YearUse,".tif")%in%paste0("./Results2/Velocity/AllModels/",RCP,"/BIOCLIM/")){
      # Step 1a. Load the ARM derived velocities        
      Vel.ARM.List <- lapply(paste0("bio",BioclimVars),
                             function(var){#(var <- names(ClimNormMn)[1])
                               rast(paste0("./Results2/Velocity/AllModels/",RCP,"/BIOCLIM/","AllModels_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
                                    lyrs = "Vel.ARM")
                             })
      
      # Step 2a. Estimate the displacement as the median velocity  
      Displ.ARM <- 10^app(log10(do.call("c",Vel.ARM.List)),
                          fun=function(i, ff) ff(i),
                          cores =3,
                          ff=function(x){median(x[x > -2 & x < 2])})
      
      # Step 1b. Load the Inter annual change derived velocities        
      Vel.IntAnn.List <- lapply(paste0("bio",BioclimVars),
                                function(var){#(var <- names(ClimNormMn)[1])
                                  rast(paste0("./Results2/Velocity/AllModels/",RCP,"/BIOCLIM/","AllModels_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
                                       lyrs = "Vel.IntAnn")
                                })
      # Step 2b. Estimate the displacement as the median velocity  
      Displ.IntAnn <- 10^app(log10(do.call("c",Vel.IntAnn.List)),
                             fun=function(i, ff) ff(i),
                             cores =3,
                             ff=function(x){median(x[x > -2 & x < 2])})
      # Step 1c. Load the Stra End Anomaly change derived velocities        
      Vel.Ann.List <- lapply(paste0("bio",BioclimVars),
                             function(var){#(var <- names(ClimNormMn)[1])
                               rast(paste0("./Results2/Velocity/AllModels/",RCP,"/BIOCLIM/","AllModels_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
                                    lyrs = "Vel.Ann")
                             })
      # Step 2b. Estimate the displacement as the median velocity  
      Displ.Ann <- 10^app(log10(do.call("c",Vel.Ann.List)),
                          fun=function(i, ff) ff(i),
                          cores =3,
                          ff=function(x){median(x[x > -2 & x < 2])})
      
      ## Estimate multivariate metrics- Divergence
      # Step 1.: Load the bearings
      BearingList <- lapply(paste0("bio",BioclimVars),
                            function(var){#(var <- names(ClimNormMn)[1])
                              rast(paste0("./Results2/Velocity/AllModels/",RCP,"/BIOCLIM/","AllModels_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
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
      OutRastMultVar <- c(Displ.ARM, Displ.IntAnn, Displ.Ann,
                          Divergence)
      names(OutRastMultVar) <- c("Displ.ARM","Displ.IntAnn","Displ.Ann","Dive")
      
      ## Save the Output
      writeRaster(OutRastMultVar,
                  filename = paste0("./Results2/Velocity/AllModels/",RCP,"/BIOCLIM/","AllModels_DispDiv_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
                  overwrite = TRUE)
      ## Clean the Memory
      rm(list = c("Displ.ARM", "Displ.IntAnn", "Displ.Ann", "Vel.ARM.List", "Vel.IntAnn.List", "Vel.Ann.List","Divergence","OutRastMultVar"));gc() 
    }
  }
}


