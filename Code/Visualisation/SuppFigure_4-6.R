rm(list=ls());gc()
require(terra)
library(tidyterra)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sp")
wrld_simpl2 <- vect(world)
wrld_simpl2 <- project(wrld_simpl2,"+proj=eck4")
rm(world);gc()

## Novelty estimated using a Mahalanobis distance 
MDAllRCP <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                   function(RCP){#(RCP <- "RCP60")
                     # Load the Treshold
                     # Load the Treshold
                     Tresh <- readRDS(paste0("./Results/Novelty/Mean_All_Models/",
                                             RCP,
                                             "/AllModels_",RCP,"_TreshSumm.rds"))
                     # Load the MD distances for each yeat                      
                     MDSum.ModlRCP <- lapply(paste0("./Results/Novelty/Mean_All_Models/",
                                                    RCP,"/AllModels_",RCP,"_",
                                                    seq(2099, 2299,by=50),"_MDminSumm.tif"),
                                             function(x){rast(x)[[1]]})
                     MDSum.ModlRCP <- do.call("c",MDSum.ModlRCP)
                     names(MDSum.ModlRCP) <- seq(2099, 2299,by=50)+1
                     return(MDSum.ModlRCP)
                   })
### Compile the Novel by when
MDAllRCP <- do.call("c",MDAllRCP)
names(MDAllRCP) <- paste0(rep(c("RCP26", "RCP45", "RCP60", "RCP85"),each=5),"-",names(MDAllRCP) )

### Plot the Climate novelty metric
pdf("./Results/PDF/Sup_Figure_4.pdf",
    width = 10, height=12)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = log10(MDAllRCP)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #color scheme 
                        na.value = NA,#Do not map NA
                        breaks = seq(-3,2,by=0.5), # Legend breaks
                        name = "Log10-Environmental\nDistance",
                        limit = c(-3,2))+ # Legend Title
  facet_wrap(~lyr, ncol=5) +
  theme(legend.position="bottom",legend.key.width = unit(4, 'cm')) +
  geom_spatvector(fill = NA,linewidth = 0.2) # Add the vector of the world
dev.off()


## Novelty estimated using a Displacement
DisplARM <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                   function(RCP){#(RCP <- "RCP26")
                     #Load the displacement for each year
                     RastTmpList <- lapply(dir(paste0("./Results/Displacement_Divergence/Mean_All_Models/",
                                                      RCP),pattern="AllModels_DispDiv_",full.names=T)[-1],
                                           function(rast_use){
                                             rast(rast_use)[[3]]
                                           })
                     RastTmp <- do.call("c",RastTmpList)
                     names(RastTmp) <- seq(2099, 2299,by=50)+1
                     return(RastTmp)
                   })
### Compile the Novel by when
DisplRCP <- do.call("c",DisplARM)
names(DisplRCP) <- paste0(rep(c("RCP26", "RCP45", "RCP60", "RCP85"),each=5),
                          "-",names(DisplRCP))

### Plot the Climate Displacement metrics 
pdf("./Results/PDF/Sup_Figure_5.pdf",
    width = 10, height=12)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = log10(DisplRCP)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = -2:2, # Legend breaks
                        labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Log10-Speed [km/yrs]", # Legend Title
                        limit = c(-2,2)) +
  facet_wrap(~lyr,ncol=5) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  theme(legend.position="bottom",legend.key.width = unit(4, 'cm'))
dev.off()

## Novelty estimated using a Displacement
Divergence <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                     function(RCP){#(RCP <- "RCP26")
                       #Load the divergecne for each year
                       RastTmpList <- lapply(dir(paste0("./Results/Displacement_Divergence/Mean_All_Models/",
                                                        RCP),pattern="AllModels_DispDiv_",full.names=T)[-1],
                                             function(rast_use){
                                               rast(rast_use)[[4]]
                                             })
                       RastTmp <- do.call("c",RastTmpList)
                       names(RastTmp) <- seq(2099, 2299,by=50)+1
                       return(RastTmp)
                     })
### Compile the Novel by when
DiveRCP <- do.call("c",Divergence)
names(DiveRCP) <- paste(rep(c("RCP26", "RCP45", "RCP60", "RCP85"),each=5),
                        "-",names(DiveRCP))
### Plot the Climate divergence metrics 
pdf("./Results/PDF/Sup_Figure_6.pdf",
    width = 10, height=12)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = DiveRCP) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = c(0,60,120,180), # Legend breaks
                        #labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Divergence angle", # Legend Title
                        limit = c(0,180)) +
  facet_wrap(~lyr,ncol = 5) +
  theme(legend.position="bottom",legend.key.width = unit(4, 'cm')) +
  geom_spatvector(fill = NA)
dev.off()
