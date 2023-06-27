rm(list=ls());gc()
require(terra)
library(tidyverse)
library(tidyterra)
require(maptools)
data(wrld_simpl)
wrld_simpl2 <- vect(wrld_simpl)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Results/Novelty/BIOCLIM2")

# Plot the Novelty for each of the evaluated time periods (2099,2149,2199,2249,2299) for all RCPs
## Novelty estimated using a Mahalanobis distance 
MDNovelYr <- lapply(seq(2099, 2299,by=50),
                    function(YearUse){#(YearUse <- seq(2099, 2299,by=50)[1])
                      MDNovelPerRCPList <- lapply(c("RCP26","RCP60","RCP85"),
                                                  function(RCP){#(RCP <- c("RCP26","RCP60","RCP85")[1])
                                                    AllRast <- lapply(dir(pattern=paste0(RCP,"_",YearUse,"_MDminSumm.tif")),
                                                                      function(x){
                                                                        resample(x = rast(x,lyrs="MD.min"),
                                                                                 y = rast(nrows=180*4, ncols=360*4))
                                                                      })
                                                    MDNovelTmp <- app(do.call("c",AllRast),median)
                                                    return(MDNovelTmp)
                                                  })#RCP <- "RCP26"
                      MDNovelPerRCP <- do.call("c",MDNovelPerRCPList)
                      names(MDNovelPerRCP) <- c("RCP26","RCP60","RCP85")
                      return(MDNovelPerRCP)
                    })
names(MDNovelYr) <- paste0("Yr",seq(2099, 2299,by=50))

## Thresholds
MDNovelTrsh <- lapply(seq(2099, 2299,by=50),
                      function(YearUse){#(YearUse <- seq(2099, 2299,by=50)[1])
                        MDNovelPerRCPList <- sapply(c("RCP26","RCP60","RCP85"),
                                                    function(RCP){#(RCP <- c("RCP26","RCP60","RCP85")[1])
                                                      AllRast <- sapply(dir(pattern=paste0(RCP,"_",YearUse,".rds")),
                                                                        function(x){
                                                                          readRDS(x)$MDSummTresh
                                                                        })
                                                      MDTreshTmp <- median(AllRast)
                                                      return(MDTreshTmp)
                                                    })#RCP <- "RCP26"
                        return(MDNovelPerRCPList)
                      })
names(MDNovelTrsh) <- paste0("Yr",seq(2099, 2299,by=50))



## Create a raster saying when a given location will become novel 
NoveltySumm <- lapply(c("RCP26","RCP60","RCP85"),
                      function(RCP){
                        RCPAllYrs <- c(MDNovelYr[[1]][[RCP]])
                        for(i in 2:length(MDNovelYr)){
                          RCPAllYrs <- c(RCPAllYrs,MDNovelYr[[i]][[RCP]])
                        }
                        names(RCPAllYrs) <- seq(2099, 2299,by=50)
                        is.RCPNvl <-  RCPAllYrs >median(MDNovelTrsh[[1]])
                        WhenNovl <- app(is.RCPNvl,
                                        function(x){
                                          if(is.na(x[[1]])){
                                            out <- NA
                                          } else{
                                            WhenNvTmp <- which(x=="TRUE")
                                            out <- ifelse(length(WhenNvTmp)==0,0,min(WhenNvTmp))
                                          }
                                          return(out)
                                        })
                        return(WhenNovl)
                      })
NoveltySumm <- do.call("c",NoveltySumm)
names(NoveltySumm) <- c("RCP26","RCP60","RCP85")

ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = as.factor(NoveltySumm)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_d(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = 0:5, # Legend breaks
                        labels = c("No Change",seq(2099, 2299,by=50)), # Legend Labels
                        name = "Novelty" # Legend Title
                        ) +
  facet_wrap(~lyr) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = "BIOCLIM2") # Fig title


