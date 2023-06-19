rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
#setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
setwd("/Volumes/Crucial X6/Data/GLOBAL/CLIMATE/FUTURE/CMIP5_Climate_copernicus/Data/Processed")
for(VarGroup in c("koeppen_geiger","BIOCLIM")){#(VarGroup <- c("koeppen_geiger","BIOCLIM")[1])
  for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[2])
    #####
    #Generate a Mean climate surface using the mean of all GCMs
    #if(length(dir(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/",VarGroup,"/")))!=320){
      TimeIn <- Sys.time()
      sfInit(parallel=TRUE, cpus=5)
      sfExport("RCP")
      sfExport("VarGroup")
      sfLibrary(terra)
      AllModSummList <- sfLapply(1980:2299,
                                 function(YearUse){#(YearUse <- c(1980:2299)[[1]])
                                   # Assess if the value exists
                                   #if(!paste0("AllModels_",RCP,"_",YearUse,".tif")%in%dir(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/",VarGroup,"/"))){
                                     # Load the Data for one year for all the relevant models
                                     AllModListPerYr <- lapply(unique(sapply(strsplit(dir(paste0("./",RCP,"/",VarGroup)),"_"),function(x){x[[1]]})),#paste0("./Data/CMIP5/Processed/",RCP,"/",VarGroup))
                                                               function(Model){#(Model <- unique(sapply(strsplit(dir(paste0("./",RCP,"/",VarGroup)),"_"),function(x){x[[1]]}))[1])
                                                                 resample(x = rast(dir(paste0("./",#"./Data/CMIP5/Processed/",
                                                                                              ifelse(YearUse<=2005,"Historical",RCP),
                                                                                              "/",VarGroup,"/"),
                                                                                       pattern = paste0(Model,"_",YearUse),
                                                                                       full.names = T)),
                                                                          y = rast(nrows=180,ncols=360))
                                                               })


                                     # Estimate the Mean for a year across all models
                                     AllModPerYrSumm <- Reduce("+", AllModListPerYr)/length(AllModListPerYr)
                                     # Save the Mean for a year summary 
                                     writeRaster(AllModPerYrSumm,
                                                 paste0("./AllModelSumm1Arcmin/",RCP,"/",VarGroup,"/AllModels_",RCP,"_",YearUse,".tif"),#"./Data/CMIP5/Processed/AllModelSumm/"
                                                 overwrite=TRUE)
                                     
                                   #}
                                    gc()
                                   return(TRUE)
                                 })
      sfStop()
      print(Sys.time()-TimeIn)
      rm(list="AllModSummList");gc()
    #}
  }  
}

VarNames <- list(c("pr","tas","tasmax","tasmin"),c("P.Seasonal","Tm.Seasonal"))

for(VarGroup in c("Monthly","Seasonal")){#(VarGroup <- c("Monthly","Seasonal")[2])
  for(VarName in VarNames[[ifelse(VarGroup=="Monthly",1,2)]]){#(VarName <- VarNames[[ifelse(VarGroup=="Monthly",1,2)]][1])
    for (RCP in c("RCP26", "RCP45", "RCP60", "RCP85")){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[2])
      #####
      #Generate a Mean climate surface using the mean of all GCMs
      #if(length(dir(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/",VarGroup,"/"), pattern = VarName))!=320){
        TimeIn <- Sys.time()
        sfInit(parallel=TRUE, cpus=5)
        sfExport("RCP")
        sfExport("VarGroup")
        sfExport("VarName")
        sfLibrary(terra)
        AllModSummList <- sfLapply(1980:2299,
                                   function(YearUse){#(YearUse <- c(1980:2299)[[1]])
                                     # Assess if the value exists
                                     #if(!paste0("AllModels_",RCP,"_",YearUse,"_",VarName,".tif")%in%dir(paste0("./Data/CMIP5/Processed/AllModelSumm/",RCP,"/",VarGroup,"/"))){
                                       # Load the Data for one year for all the relevant models
                                       AllModListPerYr <- lapply(unique(sapply(strsplit(dir(paste0("./",RCP,"/",VarGroup)),"_"),function(x){x[[1]]})),
                                                                 function(Model){#(Model <- unique(sapply(strsplit(dir(paste0("./Data/CMIP5/Processed/",RCP,"/",VarGroup)),"_"),function(x){x[[1]]}))[1])
                                                                   resample(x = rast(paste0("./",
                                                                                            ifelse(YearUse<=2005,"Historical",RCP),
                                                                                            "/",VarGroup,"/",Model,"_",YearUse,"_",VarName,".tif")),
                                                                            y = rast(nrows=180,ncols=360))
                                                                 })
                                       # Estimate the Mean for a year across all models
                                       AllModPerYrSumm <- Reduce("+", AllModListPerYr)/length(AllModListPerYr)
                                       # Save the Mean for a year summary 
                                       writeRaster(AllModPerYrSumm,
                                                   paste0("./AllModelSumm1Arcmin/",RCP,"/",VarGroup,"/AllModels_",RCP,"_",YearUse,"_",VarName,".tif"),
                                                   overwrite=TRUE)
                                       
                                     #}
                                     return(TRUE)
                                   })
        sfStop()
        print(Sys.time()-TimeIn)
      #}
        rm(list="AllModSummList");gc()
    }    
  }
}