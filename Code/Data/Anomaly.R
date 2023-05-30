rm(list=ls());gc()
require(terra)
require(rnaturalearth)
world <- vect(ne_countries(scale = "large", returnclass = "sp"))

setwd("/Users/au467796/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
ModelsAll <- c("bcc-csm1-1", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0",
               "GISS-E2-H", "GISS-E2-R", "HadGEM2-ES", "MPI-ESM-LR", "NorESM1-M")

# Create a climate normal surface for each climate model
ChelsaMask <- rast("./Data/CHELSA/Data/pr/pr_01_1979_V.2.1.tif")
ChelsaMask <- mask(ChelsaMask, world)

for (Model in c(ModelsAll)){#(Model <- ModelsAll[[1]])
  for(VarUse in c("pr","tas","tasmin","tasmax")){#(VarUse <- c("pr","tas","tasmin","tasmax")[1])
  }
}

# Load the 1960 to 2005 Historical data 
    HistFiles1980_2005 <-lapply(paste0(Model,"_",1980:2005,"_",VarUse,".tif"),
                                function(x){#(x<- paste0(Model,"_",1960:2005,"_",VarUse,".tif")[1])
                                  rast(paste0("./Data/CMIP5/Processed/Historical/Monthly/",x))
                                })
# make the CLimate normal file by adding the 2006 to 2010 Historical data (from each RCP) to HistFiles1980_2005
    ClimNorm4RCP <- lapply(c("RCP26", "RCP46", "RCP60", "RCP85"),
                                 function(RCP){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[1])
                                   RCPList <- lapply(paste0(Model,"_",2006:2010,"_",VarUse,".tif"),
                                                     function(x){#(x<- paste0(Model,"_",2006:2010,"_",VarUse,".tif")[1])
                                                     if(x%in%dir(paste0("./Data/CMIP5/Processed/",RCP,"/Monthly"))){
                                                       out <- rast(paste0("./Data/CMIP5/Processed/",RCP,"/Monthly/",x))
                                                     } else{
                                                       out <- NA
                                                     }
                                                     return(out)
                                                   })
                                   # Merge Historical and RCP lists
                                   ClimNorm4RCPTmp <- c(HistFiles1980_2005,
                                                        RCPList)
                                   # Estimate the mean of each Month for the Climate Normal
                                   ClimNorm4RCPTmp <- Reduce("+",ClimNorm4RCPTmp)/length(ClimNorm4RCPTmp) 
                                   # Return the Mean value
                                   return(ClimNorm4RCPTmp)
                               })
    names(ClimNorm4RCP) <- c("RCP26", "RCP46", "RCP60", "RCP85")
    
    Anomaly <- lapply(c("RCP26", "RCP46", "RCP60", "RCP85"),
                      function(RCP){#(RCP <- c("RCP26", "RCP46", "RCP60", "RCP85")[1])
                        for(Year in 1960:2299){#(Year <- c(1960:2299)[1])
                          x <- paste0(Model,"_",Year,"_",VarUse,".tif")
                          if(Year < 2006){
                            if(x%in%dir(paste0("./Data/CMIP5/Processed/Historical/Monthly"))){
                              RastIn <- rast(paste0("./Data/CMIP5/Processed/Historical/Monthly/",x))
                              # Estimate the anomly as Target - Clim Nomral so that Negative values mean decrese in tem, and positive means Increases
                              AnomalyTmp <- RastIn - ClimNorm4RCP[[RCP]]
                              AnomalyTmp2 <- mask(AnomalyTmp, world)
                              writeRaster(AnomalyTmp2,
                                          paste0("./Data/CMIP5/Processed/Anomalies/",Model,"/",RCP,"/",
                                                 Model,"_",Year,"_",VarUse,"_Anom.tif"),
                                          overwrite = T )
                              a <- Sys.time()
                              AnomalyTmpInterpolated <- resample(x = AnomalyTmp[[1]],
                                                                 y = ChelsaMask,
                                                                 filename = paste0(getwd(),"/Data/CMIP5/Processed/Anomalies/",Model,"/",RCP,"/",
                                                                                   Model,"_",Year,"_",VarUse,"_AnomResmp2.tif"),
                                                                 overwrite = T)
                              
                              Sys.time() - a
                              a <- Sys.time()
                              gdalUtils::gdalwarp(srcfile = paste0(getwd(),"/Data/CMIP5/Processed/Anomalies/",Model,"/",RCP,"/",
                                                                   Model,"_",Year,"_",VarUse,"_Anom.tif"),
                                                  dstfile = paste0(getwd(),"/Data/CMIP5/Processed/Anomalies/",Model,"/",RCP,"/",
                                                                   Model,"_",Year,"_",VarUse,"_AnomResmp.tif"),
                                                  tr = res(ChelsaMask),
                                                  r = "bilinear")
                              Sys.time() - a
                            } 
                          }
                        }
                      }
