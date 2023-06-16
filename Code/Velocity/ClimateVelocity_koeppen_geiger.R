rm(list=ls());gc()
require(terra)
require(snowfall)
require(maptools)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
source("./Code/Velocity/VelocityFnc.R")

# Biomes Map
BIOMES <- vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")

ModelsAll <- c("bcc-csm1-1", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0",
               "GISS-E2-H", "GISS-E2-R", "HadGEM2-ES", "MPI-ESM-LR", "NorESM1-M")


for (Model in c(ModelsAll)){#(Model <- ModelsAll[[1]])
  # Load the 1960 to 2005 Historical data 
  HistFiles1980_2005 <-lapply(paste0(Model,"_",1980:2005,"_koeppen_geigerVars.tif"),
                              function(x){#(x<- paste0(Model,"_",1960:2005,"_koeppen_geigerVars.tif")[1])
                                rast(paste0("./Data/CMIP5/Processed/Historical/koeppen_geiger/",x))
                              })
  
  # Estimate Novelty based on Future to - CLimate Normal distance
  for (RCP in c("RCP26", "RCP46", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[1])
    if(length(grep(Model,dir(paste0("./Data/CMIP5/Processed/",RCP,"/koeppen_geiger"))))!=0){
      # Load the 2006 to 2010 Historical data for an RCP to  Create a climate normal raster for each evaluated variable  
      RCPList <- lapply(paste0(Model,"_",2006:2010,"_koeppen_geigerVars.tif"),
                        function(x){#(x<- paste0(Model,"_",2006:2010,"_koeppen_geigerVars.tif")[1])
                          rast(paste0("./Data/CMIP5/Processed/",RCP,"/koeppen_geiger/",x))
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
      names(ClimNormMn) <- names(ClimNorm4RCPTmp[[1]])
      #Crop oceans
      ClimNormMn <- mask(ClimNormMn,BIOMES[BIOMES$BIOME<15,])
      # Clen the Memeory
      rm(list = ls()[!ls()%in%c("Model","RCP","koeppen_geigerVars","BIOMES","ModelsAll",
                                "HistFiles1980_2005","ClimNormMn",
                                "TempGradFnc", "SpatHetFnc","VelocityFnc","BearingFnc")]);gc()
      for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[3])
        # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
        RCPList <- lapply(paste0(Model,"_",c(2006:YearUse),"_koeppen_geigerVars.tif"),
                          function(x){#(x<- paste0(Model,"_",2006:2010,"_",VarUse,".tif")[1])
                            rast(paste0("./Data/CMIP5/Processed/",RCP,"/koeppen_geiger/",x))})
        ### Velocity per variable
        for(var in names(HistFiles1980_2005[[1]])){#(var <- names(HistFiles1980_2005[[1]])[1])
          if(!paste0(Model,"_",var,"_",RCP,"_1980to",YearUse,".tif")%in%dir("./Results/Velocity/koeppen_geiger/")){
            # Step 0: create a raster with the time series for a given variable 
            TimeSerRast <- do.call("c",
                                   lapply(c(HistFiles1980_2005,RCPList),
                                          function(x){x[[var]]}))
            TimeSerRast2 <- mask(TimeSerRast,BIOMES[BIOMES$BIOME<15,])
            # Step 1a: Estimate the Temporal gradient as the slope of the time series
            TimeHetARM <- app(TimeSerRast2,
                              fun=function(i, ff) ff(i),
                              cores = 3,
                              ff = function(x){
                                if(!is.na(x[1])){
                                  tmpDtaFrm <- data.frame(prop = x,
                                                          Time = c(1:length(x)))
                                  TimMod <- nlme::gls(prop~Time, data = tmpDtaFrm,
                                                      correlation = nlme::corARMA(p=1), 
                                                      method ="ML",
                                                      control=list(opt = "optim",
                                                                   msMaxIter = 100000))
                                  coef(summary(TimMod))["Time",c("Value")]
                                } else{
                                  NA
                                }
                              })
            
            #  Step 1b: Estimate the Temporal gradient as the Median of inter-anual changes        
            TimeHetIntAnn<- app(TimeSerRast2,
                                fun=function(i, ff) ff(i,method = "Anomaly1"),
                                cores = 3,
                                ff = TempGradFnc)
            # Step 1c: Estimate the Temporal gradient as the anomaly divided by time
            TimeHetAnom <- abs(TimeSerRast2[[dim(TimeSerRast2)[3]]]-TimeSerRast2[[1]])
            
            # Step 2.: Estimate the Spatial gradient -  Here we use the ClimNormMn raster
            SpatHetTmp <- SpatHetFnc(ClimNormMn[[var]])
            
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
            # clean memory
            rm(list = c("TimeSerRast",
                        "TimeHetARM","TimeHetIntAnn",
                        "SpatHetTmp",
                        "VelocityARM","VelocityIntAnn",
                        "BearingTmp",
                        "OutRast"))
            gc()
          }
        }
        
        ### Estimate multivariate metrics - Displacement
        # Step 1a. Load the ARM derived velocities        
        Vel.ARM.List <- lapply(names(HistFiles1980_2005[[1]]),
                               function(var){#(var <- names(HistFiles1980_2005[[1]])[1])
                                 rast(paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_1980to",YearUse,".tif"),
                                      lyrs = "Vel.ARM")
                               })
        # Step 2a. Estimate the displacement as the median velocity  
        Displ.ARM <- 10^app(log10(do.call("c",Vel.ARM.List)),
                            fun=function(i, ff) ff(i),
                            cores =3,
                            ff=function(x){median(x[x > -2 & x < 2])})
        
        # Step 1b. Load the Inter annual change derived velocities        
        Vel.IntAnn.List <- lapply(names(HistFiles1980_2005[[1]]),
                                  function(var){#(var <- names(HistFiles1980_2005[[1]])[1])
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
        BearingList <- lapply(names(HistFiles1980_2005[[1]]),
                              function(var){#(var <- names(HistFiles1980_2005[[1]])[1])
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



