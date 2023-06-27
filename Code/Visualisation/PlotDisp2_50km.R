rm(list=ls());gc()
require(terra)
library(tidyverse)
library(tidyterra)
require(maptools)
data(wrld_simpl)
wrld_simpl2 <- vect(wrld_simpl)
setwd("/Volumes/MacPro 2013 Backup/BiomeChange/")

# Plot the Novelty for each of the evaluated time periods (2099,2149,2199,2249,2299) for all RCPs
## Novelty estimated using a Mahalanobis distance 

Model <- c("BIOCLIM","BIOCLIM2","koeppen_geiger","Seasonal")[1]#(Model <- c("BIOCLIM","BIOCLIM2","koeppen_geiger","Seasonal")[1])
MDAllRCP <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                   function(RCP){#(RCP <- "RCP26")
                     # Load the Treshold
                     Tresh <- readRDS(paste0("./Results2/Novelty/AllModels_50km/",
                                             RCP,"/",Model,
                                             "/AllModels_",RCP,"_TreshSumm.rds"))
                     # Load the MD distances for each yeat                      
                     Tresh <- readRDS(paste0("./Results2/Novelty/AllModels_50km/",
                                             RCP,"/",Model,
                                             "/AllModels_",RCP,"_TreshSumm.rds"))
                     
                     MDSum.ModlRCP <- lapply(paste0("./Results2/Novelty/AllModels_50km/",
                                                    RCP,"/",Model,"/AllModels_",RCP,"_",
                                                    seq(2099, 2299,by=50),"_MDminSumm.tif"),
                                             function(x){rast(x,lyrs = 1)})
                     MDSum.ModlRCP <- do.call("c",MDSum.ModlRCP)
                     names(MDSum.ModlRCP) <- seq(2099, 2299,by=50)+1
                     # Assess if the MD corssess the treshold
                     NoAnalogue <- MDSum.ModlRCP>Tresh$MDSummTresh
                     NoAnalYr <- app(NoAnalogue,
                                     function(x){
                                       if(is.na(x[1])){
                                         NA
                                       } else{
                                         out <- seq(2099, 2299,by=50)[min(which(x==1))]
                                         ifelse(is.na(out),0,out)
                                       }
                                     })
                     return(NoAnalYr)
                   })
# Compile the Novel by when
MDAllRCP <- do.call("c",MDAllRCP)
names(MDAllRCP) <- c("RCP26", "RCP45", "RCP60", "RCP85")


ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = as.factor(MDAllRCP)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_d(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = c(0,seq(2099, 2299,by=50)), # Legend breaks
                        labels = c("No Change",seq(2099, 2299,by=50)), # Legend Labels
                        name = "Year when\nNo-Analogue?" # Legend Title
  ) +
  facet_wrap(~lyr) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = Model) # Fig title


# Proportion of the earth that will be novel by 2300
colSums(apply(values(MDAllRCP,na.rm=T),
              2,
              function(x){out <- round(table(factor(x,c(0,seq(2099, 2299,by=50))))/sum(table(x))*100,2)})[-1,])

# Barplot shown the miant of the earth to be novel 
MDAllRCPSummTbl <- lapply(1: dim(MDAllRCP)[3],
                          function(x){
                            out <- as.data.frame(table(factor(values(MDAllRCP[[x]],na.rm=T),c(0,seq(2099, 2299,by=50))))/sum(table(values(MDAllRCP[[x]],na.rm=T))))[-1,]
                            out$Freq <- round(out$Freq*100,2)
                            names(out) <- c("Year","Freq")
                            return(out)
                          })
MDAllRCPSummTbl <- data.frame(RCP = rep(names(MDAllRCP),each=dim(MDAllRCPSummTbl[[1]])[1]),
                              do.call("rbind",MDAllRCPSummTbl))
MDAllRCPSummTbl$Year <- factor(MDAllRCPSummTbl$Year,
                               levels=rev(seq(2099, 2299,by=50)))


ggplot(MDAllRCPSummTbl, aes(fill=Year, y=Freq, x=RCP)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "RdYlBu") +
  ylab("Proportion Land considered novel") +
  xlab("RCP Scenario") +
  ggtitle(Model)


