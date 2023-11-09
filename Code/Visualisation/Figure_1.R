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

# Estimate by when a given area will be novel - Novelty estimated using a Mahalanobis distance 
MDAllRCP <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                   function(RCP){#(RCP <- "RCP45")
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
                     # Assess if the MD corssess the treshold
                     NoAnalogue <- (MDSum.ModlRCP>Tresh$MDSummTresh)*1
                     NoAnalYr <- app(NoAnalogue,
                                     function(x){
                                       if(is.na(x[1])){
                                         out <- NA
                                       } else{
                                         if(sum(x)==0){
                                           out <- 0
                                         } else{
                                           out <- seq(2099, 2299,by=50)[min(which(x==1))]
                                         }
                                       }
                                       return(out)
                                     })
                     return(NoAnalYr)
                   })
# Compile the Novel by when
MDAllRCP <- do.call("c",MDAllRCP)
names(MDAllRCP) <- c("RCP26", "RCP45", "RCP60", "RCP85")
# Turn into novel area (Novel at any point)
MDAllRCP <- (MDAllRCP!=0)*1


# Estimate the mean Displacement over the evaluated time periods (2049,2099,2149,2199,2249,2299) for all RCPs - Displacement estimated using an ARM model
DisplARM <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                   function(RCP){#(RCP <- "RCP26")
                     #Load the displacement for each year
                     RastTmpList <- lapply(dir(paste0("./Results/Displacement_Divergence/Mean_All_Models/",
                                                      RCP),pattern="AllModels_DispDiv_",full.names=T)[-1],
                                           function(rast_use){
                                             rast(rast_use)[[3]]
                                           })
                     RastTmp <- do.call("c",RastTmpList)
                     names(RastTmp) <- paste(c(seq(2049, 2299,by=50)-50),"to",seq(2049, 2299,by=50))[-1]
                     return(RastTmp)})
names(DisplARM) <- c("RCP26", "RCP45", "RCP60", "RCP85")

## Make a mean summary 
DisplARMList <- lapply(DisplARM,function(x){mean(x)})
DisplARMMeanRast <- c(DisplARMList[[1]],
                      DisplARMList[[2]],
                      DisplARMList[[3]],
                      DisplARMList[[4]])
names(DisplARMMeanRast) <- c("RCP26", "RCP45", "RCP60", "RCP85")
# Turn into novel area (Displacment over 1)
DisplARMMeanRast <- (DisplARMMeanRast > 1)*1


# Estimate the mean divergence accorss evaluated time periods (2049,2099,2149,2199,2249,2299) for all RCPs
# Displacement estimated using an ARM model

Divergence <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                     function(RCP){#(RCP <- "RCP26")
                       #Load the displacement for each year
                       RastTmpList <- lapply(dir(paste0("./Results/Displacement_Divergence/Mean_All_Models/",
                                                        RCP),pattern="AllModels_DispDiv_",full.names=T)[-1],
                                             function(rast_use){
                                               rast(rast_use)[[4]]
                                             })
                       RastTmp <- do.call("c",RastTmpList)
                       names(RastTmp) <- paste(c(seq(2049, 2299,by=50)-50),"to",seq(2049, 2299,by=50))[-1]
                       return(RastTmp)})
names(Divergence) <- c("RCP26", "RCP45", "RCP60", "RCP85")

## Make a mean summary 
DivergenceSumm <- c(mean(Divergence[[1]]),
                    mean(Divergence[[2]]),
                    mean(Divergence[[3]]),
                    mean(Divergence[[4]]))
names(DivergenceSumm) <- names(Divergence)
# Turn into novel area (ortogonal Divergence)
DivergenceSumm <- (DivergenceSumm>60) * (DivergenceSumm<120)

# Summary of all novelty
ContTbl <- rbind(c(0,0,0),#31a354
                 c(1,0,0),#756bb1
                 c(0,1,0),#bcbddc
                 c(0,0,1),#efedf5
                 c(1,1,0), #ffffb2
                 c(0,1,1),#ffeda0
                 c(1,0,1),#feb24c
                 c(1,1,1))#e41a1c
#ContTbl <- expand.grid(c(0,1),c(0,1),c(0,1))

# The distance raster might be in a diff extent that the velocity raster so to ensure they match
MDAllRCP2 <- lapply(MDAllRCP,function(x){
  resample(x,DisplARMMeanRast[[1]])
  })
names(MDAllRCP2) <- c("RCP26", "RCP45", "RCP60", "RCP85")
                   
# MErge novelty maps
NovelCriteriaList <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                            function(x){#(x<-"RCP26")
                              c(MDAllRCP2[[x]],
                                DisplARMMeanRast[[x]],
                                DivergenceSumm[[x]])
                            })
NovelCriteriaList[[1]]
# Assess by which measure a given are is novel.
PlotTestList <- lapply(NovelCriteriaList,
                       function(Rast){
                         app(Rast,
                             function(x){
                               out <- which(colSums(apply(ContTbl,1,function(y){y==x}))==3)
                               out <- ifelse(length(out)!=0,out,NA)
                               return(out)
                             })
                       })

PlotTest <- do.call("c",PlotTestList)
names(PlotTest) <- c("RCP26", "RCP45", "RCP60", "RCP85")
PlotTest <- as.factor(PlotTest)

round(apply(values(PlotTest,na.rm=T),2,function(x){table(factor(x,levels=1:8))/length(x)})*100,3)

pdf("./Results/PDF/Figure_1.pdf",
    width = 5, height=7)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = PlotTest) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_manual(values = c("#31a354","#54278f","#756bb1","#cbc9e2","#ffffb2","#ffeda0","#feb24c","#e41a1c"),
                    na.value = NA,#Do not map NA
                    breaks = 1:8, # Legend breaks
                    labels = c("No Novelty",
                               "NovClim",
                               "FastDis",
                               "OrthDiv",
                               "NovClim/FastDis",
                               "FastDis/OrthDiv",
                               "NovClim/OrthDiv",
                               "NovClim/FastDis/OrthDiv"), # Legend Labels
                    name = "Source of Novelty" # Legend Title
  ) +
  facet_wrap(~lyr, ncol=1) +
  geom_spatvector(fill = NA,linewidth = 0.2) # Add the vector of the world
  
dev.off()
