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
  # Estimate Novelty based on Future to - Climate Normal distance
  for (RCP in c("RCP26", "RCP46", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[1])
    for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[1])
      if(length(grep(Model,dir(paste0("./Data/CMIP5/Processed/",RCP,"/koeppen_geiger"))))!=0){
        #Create a climate normal raster for each evaluated variable  based on the data for the 30 years before the last start period 
        ClimNorm4RCPTmp <- lapply(paste0(Model,"_",c(YearUse-80):c(YearUse-50),"_koeppen_geigerVars.tif"),
                                  function(x){#(x<- paste0(Model,"_",2006:2010,"_koeppen_geigerVars.tif")[1])
                                    rast(paste0("./Data/CMIP5/Processed/",RCP,"/koeppen_geiger/",x))
                                  })
        # Estimate the mean of each band for the Climate Normal
        ClimNormMn <- do.call("c",
                              lapply(names(ClimNorm4RCPTmp[[1]]),
                                     function(VarUse){#(VarUse <- names(ClimNorm4RCPTmp[[1]])[1])
                                       app(do.call("c",
                                                   lapply(ClimNorm4RCPTmp,
                                                          function(x){x[[VarUse]]})),mean)
                                     }))
        
        names(ClimNormMn) <- c("T_w.m", "T_c.m","T_avg","P_tot","P_wint","P_summ","P_d.m","P_d.m.summ","P_d.m.wint","P_w.m","P_w.m.summ","P_w.m.wint","T_4th_w.m")
        #Crop oceans
        ClimNormMn <- mask(ClimNormMn,BIOMES[BIOMES$BIOME<15,])
        # Clean the Memory
        rm(list = ls()[!ls()%in%c("Model","RCP","koeppen_geigerVars","BIOMES","ModelsAll",
                                  "HistFiles1980_2005","ClimNormMn","YearUse",
                                  "TempGradFnc", "SpatHetFnc","VelocityFnc","BearingFnc")]);gc()
        # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
        RCPList <- lapply(paste0(Model,"_",c(c(YearUse-50):YearUse),"_koeppen_geigerVars.tif"),
                          function(x){#(x<- paste0(Model,"_",2006:2010,"_",VarUse,".tif")[1])
                            rast(paste0("./Data/CMIP5/Processed/",RCP,"/koeppen_geiger/",x))})
        ### Velocity per variable
        for(var in names(ClimNormMn)){#(var <- names(ClimNormMn)[1])
          if(!paste0(Model,"_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif")%in%dir("./Results/Velocity/koeppen_geiger/")){
            # Step 0: create a raster with the time series for a given variable 
            TimeSerRast <- do.call("c",
                                   lapply(RCPList,
                                          function(x){x[[var]]}))
            # Step 1a: Estimate the Temporal gradient as the slope of the time series
            TimeSerRast2 <- mask(TimeSerRast,BIOMES[BIOMES$BIOME<15,])
            TimeHetARM <- app(TimeSerRast2,
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
                                }
                                if(length(unique(x))==1){
                                  0
                                }
                                else{
                                  NA
                                }
                              })
            #  Step 1a: Estimate the Temporal gradient as the Median of inter-annual changes        
            TimeHetIntAnn<- app(TimeSerRast,
                                fun=function(i, ff) ff(i,method = "Anomaly1"),
                                cores = 3,
                                ff = TempGradFnc)
            # Step 1.: Estimate the Spatial gradient -  Here we use the ClimNormMn raster
            SpatHetTmp <- SpatHetFnc(ClimNormMn[[var]])
            # Step 3.: Estimate the Velocity as the ratio between spatial and temporal gradients 
            # estimate the Velocity using a slope of the time series derived Temporal gradient
            VelocityARM <- VelocityFnc(TimeHetARM,
                                       SpatHetTmp)
            # estimate the Velocity using a sMedian of inter-anual changes derived Temporal gradient
            VelocityIntAnn <- VelocityFnc(TimeHetIntAnn,
                                          SpatHetTmp)
            # Step 4.: Estimate the Bearing using the degrees from north of the vector-sum used to estimate the Spatial gradient
            BearingTmp <- BearingFnc(c(ClimNormMn[[var]],
                                       RCPList[[length(RCPList)]][[var]]))
            # Step 5.: Create a summary for a variable
            OutRast <- c(TimeHetARM,TimeHetIntAnn,SpatHetTmp,BearingTmp,VelocityARM,VelocityIntAnn)
            names(OutRast) <- c("TimeH.ARM", "TimeH.IntAnn", "SpatHet", "Bearing", "Vel.ARM", "Vel.IntAnn")
            #Save the Output
            writeRaster(OutRast,
                        filename = paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
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
        Vel.ARM.List <- lapply(names(ClimNormMn),
                               function(var){#(var <- names(ClimNormMn)[1])
                                 rast(paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
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
                                    rast(paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
                                         lyrs = "Vel.IntAnn")
                                  })
        # Step 2b. Estimate the displacement as the median velocity  
        Displ.IntAnn <- 10^app(log10(do.call("c",Vel.IntAnn.List)),
                               fun=function(i, ff) ff(i),
                               cores =3,
                               ff=function(x){median(x[x > -2 & x < 2])})
        ## Estimate multivariate metrics- Divergence
        # Step 1.: Load the bearings
        BearingList <- lapply(names(ClimNormMn),
                              function(var){#(var <- names(ClimNormMn)[1])
                                rast(paste0("./Results/Velocity/koeppen_geiger/",Model,"_",var,"_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
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
        OutRastMultVar <- c(Displ.ARM, Displ.IntAnn,Divergence)
        names(OutRastMultVar) <- c("Displ.ARM","Displ.IntAnn","Dive")
        
        ## Save the Output
        writeRaster(OutRastMultVar,
                    filename = paste0("./Results/Displacement_Divergence/koeppen_geiger/",Model,"_DispDiv_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
                    overwrite = TRUE)
        ## Clean the Memory
        rm(list = ls()[!ls()%in%c("Model","RCP","koeppen_geigerVars","BIOMES","ModelsAll",
                                  "HistFiles1980_2005","ClimNormMn",
                                  "TempGradFnc", "SpatHetFnc","VelocityFnc","BearingFnc")]);gc()
      }
    }
  }
}



