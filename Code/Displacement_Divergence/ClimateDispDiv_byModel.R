rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
setwd("~/Documents/GitHub/BiomeChange")
source("./Code/Displacement_Divergence/VelocityFnc.R")

# Define the "BIOME" for each location
BIOMES <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_100km.tif")

#####
# Estimate Novelty based on Climate Velocity metrics
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")[3:4]){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[3])
  for(ModelUse in c("bcc-csm1-1", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0",
                    "GISS-E2-H", "GISS-E2-R", "HadGEM2-ES","MPI-ESM-LR", "NorESM1-M")){#(ModelUse <- "CCSM4")
    for(YearUse in seq(2050,2300,by=50)){#(YearUse <- seq(2050,2300,by=50)[2])
      # Ensure there is future data and that the Disimilarity has not been estimated
      if(#!paste0("Displace_Diverg_",ModelUse,"_",RCP,"_",YearUse,".tif")%in% dir(paste0("/Users/au467796/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Results2/Velocity/PerModel_Novelty/",RCP,"/Seasonal/")) &
         length(dir(paste0("./",RCP,"/Seasonal"),pattern = paste0(ModelUse,"_",YearUse-15)))!=0)
        {
        TimeIn2 <- Sys.time()
        # Load the tart-Period Climatic conditions 
        if(YearUse==2050){
          StartPeriod <- c(rast(dir("./Historical/Seasonal",pattern=ModelUse, full.names = T)[1]),
                           rast(dir("./Historical/Seasonal",pattern=ModelUse, full.names = T)[2]))
        } else {
          StartPeriod <- c(rast(dir(paste0("./",RCP,"/Seasonal"),pattern = paste0(ModelUse,"_",YearUse-65),full.names = T)[1]),
                           rast(dir(paste0("./",RCP,"/Seasonal"),pattern = paste0(ModelUse,"_",YearUse-65),full.names = T)[2]))
        }
        # Project to eacal area
        StartPeriodPrj <- project(StartPeriod,"+proj=eck4",method="near",threads=50)
        # resample to 50 x 50 km
        StartPeriod_50x50 <- lapply(1:dim(StartPeriodPrj)[3],
                                    function(x){
                                      round(resample(StartPeriodPrj[[x]],eck4Rast,threads=50),ifelse(x>=9,0,1))})
        rm(list=c("StartPeriod","StartPeriodPrj"));gc()
        # Step 1.: Estimate the Spatial gradient -  Here we use the StartPeriod_50x50 raster
        SpatHetTmp <- lapply(StartPeriod_50x50,
                             function(x){
                               SpatHetFnc(RastIn = x,
                                          Dist = res(x)[1]/1000)
                             })
        
        # Load the End-Period CLimatic conditions 
        EndPeriod <- c(rast(dir(paste0("./",RCP,"/Seasonal"),pattern = paste0(ModelUse,"_",YearUse-15),full.names = T)[1]),
                       rast(dir(paste0("./",RCP,"/Seasonal"),pattern = paste0(ModelUse,"_",YearUse-15),full.names = T)[2]))
        # Project to eacal area
        EndPeriodPrj <- project(EndPeriod,"+proj=eck4",method="near",threads=50)
        # resample to 50 x 50 km
        EndPeriod_50x50 <- lapply(1:dim(EndPeriodPrj)[3],
                                  function(x){
                                    round(resample(EndPeriodPrj[[x]],eck4Rast,threads=50),ifelse(x>=9,0,1))})
        rm(list=c("EndPeriod","EndPeriodPrj"));gc()
        #  Step 2: Estimate the Temporal gradient as the Dif between Start and end conditions        
        TimeHetAnomaly <- lapply(1:length(EndPeriod_50x50),
                                 function(x){
                                   (EndPeriod_50x50[[x]]-StartPeriod_50x50[[x]])/50
                                 })
        # Step 3: Estimate the Velocity as the ratio between spatial and temporal gradients 
        VelocityAnomaly <- lapply(1:length(EndPeriod_50x50),
                                  function(x){
                                    VelocityFnc(TimeHetAnomaly[[x]],
                                                SpatHetTmp[[x]])
                                  })
        # Step 4: Estimate the Bearing using the degrees from north of the vector-sum used to estimate the Spatial gradient
        BearingTmp <- lapply(1:length(EndPeriod_50x50),
                             function(x){
                               BearingFnc(c(StartPeriod_50x50[[x]],
                                            EndPeriod_50x50[[x]]))
                             })
        # Step 5.: Create a summary for each variable and save it
        OutRast <- c(do.call("c",VelocityAnomaly),
                     do.call("c",BearingTmp))
        names(OutRast) <- paste0(rep(c("Vel.","Bear."),each=8),
                                 rep(paste0(rep(c("P.","T."),each=4),
                                            rep(c("MAP", "JJA", "SON", "DJF"),2)),2))
        #Save the Output
        writeRaster(OutRast,
                    filename = paste0("/Users/alejandroordonez/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Results2/Velocity/PerModel_Novelty/",
                                      RCP,"/Seasonal/Velocity_Bering_",ModelUse,"_",RCP,"_",YearUse,".tif"),
                    overwrite = TRUE)
        rm(list=c("EndPeriod_50x50","StartPeriod_50x50", "TimeHetAnomaly","SpatHetTmp","OutRast"));gc()
        print(Sys.time()-TimeIn2)
        # Step 6. Estimate the displacement as the median velocity  
        Displ.Anomaly <- 10^app(log10(do.call("c",VelocityAnomaly)),
                                fun=function(i, ff) ff(i),
                                cores =3,
                                ff=function(x){median(x[x > -2 & x < 2])})
        # Step 7. Estimate the Divergence as the median difference between bearings -  median angle between bearings
        Divergence <- app(do.call("c",BearingTmp),
                          function(x){
                            if(is.na(x[1])){out<-NA
                            } else{x <- x[x!=361]
                            #x <- round(x,1) 
                            y <- dist(na.omit(x))
                            y[y>180] <- 360 - y[y>180]
                            out <- median(y) }
                            return((out))
                          })
        # Step 8. Save the Multivariate Output
        OutRastMultVar <- c(Displ.Anomaly,
                            Divergence)
        names(OutRastMultVar) <- c("Displ.Ann","Dive")
        #Save the Output
        writeRaster(OutRastMultVar,
                    filename = paste0("/Users/alejandroordonez/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Results2/Velocity/PerModel_Novelty/",
                                      RCP,"/Seasonal/Displace_Diverg_",ModelUse,"_",RCP,"_",YearUse,".tif"),
                    overwrite = TRUE)
        rm(list=c("VelocityAnomaly","BearingTmp","OutRastMultVar","Displ.Anomaly", "Divergence"));gc()
        print(Sys.time()-TimeIn2)
      }
      print(YearUse)
    }
    print(ModelUse)
  }
  print(RCP)
}