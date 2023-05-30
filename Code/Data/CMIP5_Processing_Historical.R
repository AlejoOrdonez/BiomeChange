rm(list=ls());gc()

# Load the required packages
require(terra)
require(ncdf4)
require(ncdf4.helpers)
require(dismo) # To estimate Bioclimatic vars
require(ClimClass) # To estimate koeppen_geiger vars and classification
#WD.Raw <- "/Volumes/Seagate Expansion/Data/GLOBAL/CLIMATE/FUTURE/CMIP5_2300_Copernicus/Raw"
WD.Raw <-       "~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Data/CMIP5/Raw"
WD.Processed <- "~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Data/CMIP5/Processed"
# Historical data

for(Model in dir(WD.Raw)){#(Model <- dir(WD.Raw)[1])#"CCSM4"
  setwd(paste0(WD.Raw,"/",Model,"/Historical"))
  # Set the target Year
  for(TargYr in 1960:2005){#(TargYr <- 1960)
    if(!paste0(Model,"_",TargYr,"_koeppen_geigerVars.tif")%in%dir(paste0(WD.Processed,"/Historical/koeppen_geiger/"))){
      a <-Sys.time()
      
      # Get the Time Series for each model -  used to define the "Layers" to use [those representing the 12 months for the selected year]
      TimSerList <- lapply(dir(pattern="pr"),
                           function(i){
                             f <- nc_open(i)
                             ts <- nc.get.time.series(f)              
                             nc_close(f)
                             return(ts)
                           })
      ts <- do.call("c",TimSerList)
      ValsYrUse <- which(format(ts,format="%Y") == TargYr) # Layers representing the 12 months for the selected year
      
      # Load climatology for a target Year - Precipitation Data
      PList <- lapply(dir(pattern="pr"),
                      function(i){
                        rast(i)
                      })
      P <- do.call("c",PList)[[ValsYrUse]]*(86400*30) # Original P in kg m-2 s-1. Multiplying it by 86400 to convert into mm/day, and then multiplying by 30 makes it a monthly total.
      crs(P) <- "epsg:4326" # Ensure data is Projected
      rm(PList);gc()
      terra::rotate(x = P,
                    filename = paste0(WD.Processed,"/Historical/Monthly/",Model,"_",TargYr,"_pr.tif"),
                    overwrite = TRUE)
      
      # Load climatology for a target Year - Minimum Temperature Data
      TnList <- lapply(dir(pattern="tasmin"),
                       function(i){
                         rast(i)
                       })
      Tn <- do.call("c",TnList)[[ValsYrUse]] - 273.15 # Original Temp in K. Subtracting 273.15 to convert Temperatures to C
      crs(Tn) <- "epsg:4326" # Ensure data is Projected
      rm(TnList);gc()
      terra::rotate(x = Tn,
                    filename = paste0(WD.Processed,"/Historical/Monthly/",Model,"_",TargYr,"_tasmin.tif"),
                    overwrite = TRUE)
      
      # Load climatology for a target Year - Maximum Temperature Data
      TxList <- lapply(dir(pattern="tasmax"),
                       function(i){
                         rast(i)
                       })
      Tx <- do.call("c",TxList)[[ValsYrUse]] - 273.15 # Original Temp in K. Subtracting 273.15 to convert Temperatures to C
      crs(Tx) <- "epsg:4326" # Ensure data is Projected
      rm(TxList);gc()
      terra::rotate(x = Tx,
                    filename = paste0(WD.Processed,"/Historical/Monthly/",Model,"_",TargYr,"_tasmax.tif"),
                    overwrite = TRUE)
      
      # Load climatology for a target Year - Mean Temperature Data
      TmList <- lapply(dir(pattern="tas_"),
                       function(i){
                         rast(i)
                       })
      Tm <- do.call("c",TmList)[[ValsYrUse]] - 273.15 # Original Temp in K. Subtracting 273.15 to convert Temperatures to C
      crs(Tm) <- "epsg:4326" # Ensure data is Projected
      rm(TmList);gc()
      terra::rotate(x = Tm,
                    filename = paste0(WD.Processed,"/Historical/Monthly/",Model,"_",TargYr,"_tas.tif"),
                    overwrite = TRUE)
      
      # Estimate BIOCLIM - Using biovars from the Dismo package 
      BICLIM.Vars <- biovars(prec = stack(P),
                             tmin = stack(Tn),
                             tmax = stack(Tx))
      BIOCLIM.Vars2 <- rast(BICLIM.Vars)
      terra::rotate(x = BIOCLIM.Vars2,
                    filename = paste0(WD.Processed,"/Historical/BIOCLIM/",Model,"_",TargYr,"_BIOCLIM.tif"),
                    overwrite = TRUE)
      
      # Estimate total Precipitation By Season
      P.Seasonal <- c(app(P, function(x){sum(x[c(3,4,5)])}),
                      app(P, function(x){sum(x[c(6,7,8)])}),
                      app(P, function(x){sum(x[c(9,10,11)])}),
                      app(P, function(x){sum(x[c(12,1,2)])}))
      names(P.Seasonal) <- c("MAP","JJA","SON","DJF")
      # Rotate and save the total Precipitation By Season
      terra::rotate(x = P.Seasonal,
                    filename = paste0(WD.Processed,"/Historical/Seasonal/",Model,"_",TargYr,"_P.Seasonal.tif"),
                    overwrite = TRUE)
      
      # Estimate Mean Temperature By Season
      Tm.Seasonal <- c(app(Tm, function(x){mean(x[c(3,4,5)])}),
                       app(Tm, function(x){mean(x[c(6,7,8)])}),
                       app(Tm, function(x){mean(x[c(9,10,11)])}),
                       app(Tm, function(x){mean(x[c(12,1,2)])}))
      names(Tm.Seasonal) <- c("MAP","JJA","SON","DJF")
      # Rotate and save the Mean Temperature By Season
      terra::rotate(x = Tm.Seasonal,
                    filename = paste0(WD.Processed,"/Historical/Seasonal/",Model,"_",TargYr,"_Tm.Seasonal.tif"),
                    overwrite = TRUE)
      
      # Estimate koeppen_geiger vars
      AllVas <- c(P, Tn, Tx , Tm) # Merge all Vars in a single SpatRaster
      # Function to estimate koeppen_geiger vars
      koeppen_geigerFnc <- function(i) {
        Tabl <- data.frame(month = 1:12,
                           P = i[c(1:12)],
                           Tn = i[c(13:24)],
                           Tx = i[c(25:36)],
                           Tm = i[c(37:48)])
        as.numeric(ClimClass::koeppen_geiger(Tabl,
                                             A_B_C_special_sub.classes=TRUE)[1,-14])
      }
      # Estimate koeppen_geiger vars using a parallel approach 
      koeppen_geigerVars <- app(AllVas, # The Merged SpatRaster
                                fun=function(i, ff) ff(i), # How to specify the function to run in parallel
                                cores = 10, # Number of cores to divide the task
                                ff = koeppen_geigerFnc) # Function to Use
      names(koeppen_geigerVars) <- c("T_w.m", "T_c.m", "T_avg", "P_tot", "P_wint", "P_summ", "P_d.m", "P_d.m.summ", "P_d.m.wint", "P_w.m", "P_w.m.summ", "P_w.m.wint", "T_4th_w.m")
      
      terra::rotate(x = koeppen_geigerVars,
                    filename = paste0(WD.Processed,"/Historical/koeppen_geiger/",Model,"_",TargYr,"_koeppen_geigerVars.tif"),
                    overwrite = TRUE)
      
      a - Sys.time()
      # Clean memory
      rm(list=ls()[!ls()%in%c("WD.Raw","WD.Processed","Model","TargYr")]);gc()  
    }     
  }
}



