rm(list=ls());gc()
require(terra)
library(tidyverse)
library(tidyterra)
require(maptools)
data(wrld_simpl)
wrld_simpl2 <- vect(wrld_simpl)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Results/Displacement_Divergence/Seasonal")
#Model <- c("bcc-csm1-1","CanESM2","CCSM4","CESM1-CAM5","GISS-E2-H","GISS-E2-R","HadGEM2-ES","MPI-ESM-LR")[2]

# Plot the median 1980-to-2100 vs 1980-to-2300 displacement for All evaluated RCPs.
## Displacement based on temporal anomalies using the median inter-annual changes.
MedDispYr <- lapply(seq(2099, 2299,by=50),
                    function(YearUse){#YearUse <- seq(2099, 2299,by=50)[1]
                      MedDispPerRCPList <- lapply(c("RCP26","RCP60","RCP85"),
                                                  function(RCP){#RCP <- c("RCP26","RCP60","RCP85")[1]
                                                    AllRast <- lapply(dir(pattern=paste0("DispDivAll_",RCP,"_1980to",YearUse)),
                                                                      function(x){
                                                                        resample(x = rast(x,lyrs="Displ.IntAnn"),
                                                                                 y = rast(nrows=180*4, ncols=360*4))
                                                                      })
                                                    MedDispl <- app(do.call("c",AllRast),median)
                                                    return(MedDispl)
                                                  })#RCP <- "RCP26"
                      MedDispPerRCP <- do.call("c",MedDispPerRCPList)
                      names(MedDispPerRCP) <- c("RCP26","RCP60","RCP85")
                      return(MedDispPerRCP)
                    })


# Plot the Velocities 2100
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = log10(MedDispYr[[1]])) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = -2:2, # Legend breaks
                        labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Speed\n[km*100yrs]", # Legend Title
                        limit = c(-2,2)) +
  facet_wrap(~lyr) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = paste0("Displacement", "1980 to ", seq(2099, 2299,by=50)[1])) # Fig title



# Plot the Velocity for each of the 5 time points under RCP8.5 
MedDisp_YrRCP <- do.call("c",lapply(MedDispYr,function(x){x[["RCP85"]]}))
names(MedDisp_YrRCP) <- seq(2099, 2299,by=50)

ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = log10(MedDisp_YrRCP)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = -2:2, # Legend breaks
                        labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Speed\n[km*100yrs]", # Legend Title
                        limit = c(-2,2)) +
  facet_wrap(~lyr) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = "Displacement RCP85") # Fig title


