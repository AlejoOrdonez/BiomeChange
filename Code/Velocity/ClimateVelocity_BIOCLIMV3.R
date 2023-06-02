rm(list=ls());gc()
require(terra)
require(snowfall)
require(maptools)
require(analogue)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
source("./Code/Velocity/VelocityFnc.R")
# BioclimVars to use
BioclimVars <- c(8, # Mean Temperature of Wettest Quarter
                 9, # Mean Temperature of Driest Quarter
                 10,# Mean Temperature of Warmest Quarter
                 11,# Mean Temperature of Coldest Quarter
                 16,# Precipitation of Wettest Quarter
                 17,# Precipitation of Driest Quarter
                 18,# Precipitation of Warmest Quarter
                 19)# Precipitation of Coldest Quarter

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
      
      # Estimate Novelty by comparing future climate to the Climate-Normal
      # Estimate the distance to each period
      rm(list = ls()[!ls()%in%c("Model","RCP","BioclimVars","BIOMES","ModelsAll",
                                "HistFiles1980_2005","ClimNormMn",
                                "TempGradFnc", "SpatHetFnc","VelocityFnc","BearingFnc")]);gc()
      for(YearUse in seq(2099,2299,by=50)){#(YearUse <- seq(2099,2299,by=50)[2])
        # Load the Future points (and 10 yr period of conditions up to the point of interest) data for an RCP
        RCPList <- lapply(paste0(Model,"_",c(2006:YearUse),"_BIOCLIM.tif"),
                          function(x){#(x<- paste0(Model,"_",2006:2010,"_",VarUse,".tif")[1])
                            rast(paste0("./Data/CMIP5/Processed/",RCP,"/BIOCLIM/",x))[[BioclimVars]]})
        
## Velocity per variable
        VelPerVarList <- lapply(names(HistFiles1980_2005[[1]]),
                                function(var){
                                   # Step 0: create a raster with the time series for a given variable 
                                   TimeSerRast <- do.call("c",
                                                          lapply(c(HistFiles1980_2005,RCPList),
                                                                 function(x){x[[var]]}))
                                   ### Step 1a: Estimate the Temporal gradient as the slope of the time series
                                   TimeHetARM <- app(TimeSerRast,
                                                     fun=function(i, ff) ff(i),
                                                     cores = 3,
                                                     ff = TempGradFnc)
                                   #  Step 1a: Estimate the Temporal gradient as the Median of inter-anual changes        
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
                                   OutList <- list(TimeHet = list(ARM = TimeHetARM, 
                                                                  IntAnn = TimeHetIntAnn),
                                                   SpatHet = SpatHetTmp,
                                                   Bearing = BearingTmp,
                                                   Velocity = list(ARM = VelocityARM,
                                                                   IntAnn = VelocityIntAnn))
                                   return(OutList)
                              })

        
                
### Step 4.: Create a final summary all Values
        Out.List <- list(RCP = RCP,
                         PerVar = VelPerVarList,
                         Summary = )
        #Save the Output
        saveRDS(VelPerVarList,
                paste0("./Results/Novelty/VelTmp",Model,"_",RCP,"_",YearUse,".rds"))
        rm(list = ls()[!ls()%in%c("Model","RCP","BioclimVars","BIOMES","ModelsAll","HistFiles1980_2005","MDtresh","SEDtresh",
                                  "ClimNormMn","ClimNormSD","YearUse")]);gc()
      }
    }
  }
}



