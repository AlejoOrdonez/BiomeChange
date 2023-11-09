rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
setwd("~/Documents/GitHub/BiomeChange")
source("./Code/Displacement_Divergence/VelocityFnc.R")

# Define the "BIOME" for each location
BIOMES <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_100km.tif")

# Estimate Novelty based on Future to - Climate Normal distance
for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[2])
  for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[1])
    ### Velocity per variable
    #  Step 1: Estimate the Temporal gradient as the Dif between Start and end conditions
    # Load the 1980 to 2010 Historical data for an RCP to Create a climate normal raster for each evaluated variable  
    if(YearUse==2099){
      ClimNormMn <- rast(paste0("./Data/Raw/Mean_All_Models/",RCP,"/ClimNormMn_1980-2010_",RCP,".tif"))
    } else {
      ClimNormMn <- rast(paste0("./Data/Raw/Mean_All_Models/",RCP,"/RCPFull_",c(YearUse-50),"_",RCP,".tif"))
    }
    #Create a Future climate raster for each evaluated variable  based on a 30 years period centered at the end of the period  
    ClimFut <- rast(paste0("./Data/Raw/Mean_All_Models/",RCP,"/RCPFull_",YearUse,"_",RCP,".tif"))
    TimeHetAnn <- (ClimFut - ClimNormMn)/50
    # Step 2.: Estimate the Spatial gradient -  Here we use the ClimNormMn raster
    SpatHetTmp <- lapply(1:dim(ClimNormMn)[3],
                         function(x){#(x<-1)
                           SpatHetFnc(RastIn = ClimNormMn[[x]],
                                      Dist = res(ClimNormMn)[1]/1000)
                         })
    SpatHetTmp <- do.call("c",SpatHetTmp)
    # Step 3: Estimate the Velocity as the ratio between spatial and temporal gradients 
    VelocityAnn <- lapply(1:dim(SpatHetTmp)[3],
                          function(x){
                            VelocityFnc(TimeHetAnn[[x]],
                                        SpatHetTmp[[x]])
                          })
    VelocityAnn <- do.call("c",VelocityAnn)
    # Step 4.: Estimate the Bearing using the degrees from north of the vector-sum used to estimate the Spatial gradient
    BearingTmp <- lapply(1:dim(SpatHetTmp)[3],
                         function(x){
                           BearingFnc(c(ClimNormMn[[x]],
                                        ClimFut[[x]]))
                         })
    BearingTmp <- do.call("c",BearingTmp)
    rm()
    rm(list = c("ClimNormMn","ClimFut","TimeHetAnn","SpatHetTmp"));gc()
    ### Estimate multivariate metrics- displacement/Divergence
    # Step 5. Estimate the displacement as the median velocity  
    Displ.Ann <- 10^app(VelocityAnn,
                        fun=function(i, ff) ff(i),
                        cores =3,
                        ff=function(x){median(x[x > -2 & x < 2])})
    # Step 6 Estimate the Divergence as the median difference between bearings -  median angle between bearings
    Divergence <- app(BearingTmp,
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
    OutRastMultVar <- c(Displ.Ann,
                        Divergence)
    names(OutRastMultVar) <- c("Displ.Ann","Dive")
    ## Save the Output
    writeRaster(OutRastMultVar,
                filename = paste0("./Results/Displacement_Divergence/Mean_All_Models/",RCP,"/AllModels_DispDiv_",RCP,"_",c(YearUse-50),"to",YearUse,".tif"),
                overwrite = TRUE)
    rm(list = c("Divergence","Displ.Ann","OutRastMultVar"));gc()
  }
}
