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

# Estimate by when a given area will be novel

## Novelty estimated using a Mahalanobis distance 
MDAllRCP <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                   function(RCP){#(RCP <- "RCP60")
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
### Compile the Novel by when
MDAllRCP <- do.call("c",MDAllRCP)
names(MDAllRCP) <- c("RCP26", "RCP45", "RCP60", "RCP85")


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
                     # Assess if the MD corssess the treshold
                     NoAnalogue <- RastTmp>1
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
### Compile the Novel by when
DisplRCP <- do.call("c",DisplARM)
names(DisplRCP) <- c("RCP26", "RCP45", "RCP60", "RCP85")

## Novelty estimated using a Displacement
Divergence <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                     function(RCP){#(RCP <- "RCP26")
                       #Load the divergence for each year
                       RastTmpList <- lapply(dir(paste0("./Results/Displacement_Divergence/Mean_All_Models/",
                                                        RCP),pattern="AllModels_DispDiv_",full.names=T)[-1],
                                             function(rast_use){
                                               rast(rast_use)[[4]]
                                             })
                       RastTmp <- do.call("c",RastTmpList)
                       names(RastTmp) <- seq(2099, 2299,by=50)+1
                       # Assess if the MD corssess the treshold
                       NoAnalogue <- ( RastTmp > 60) * (RastTmp < 120)
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
### Compile the Novel by when
DiveRCP <- do.call("c",Divergence)
names(DiveRCP) <- c("RCP26", "RCP45", "RCP60", "RCP85")


##### Figure 2
# Estimate the percentage change of the area of all areas for each evaluate year
## Climate Novelty
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
## Displacement Novelty
DisplRCPSummTbl <- lapply(1: dim(DisplRCP)[3],
                          function(x){
                            out <- as.data.frame(table(factor(values(DisplRCP[[x]],na.rm=T),c(0,seq(2099, 2299,by=50))))/sum(table(values(DisplRCP[[x]],na.rm=T))))[-1,]
                            out$Freq <- round(out$Freq*100,2)
                            names(out) <- c("Year","Freq")
                            return(out)
                          })
DisplRCPSummTbl <- data.frame(RCP = rep(names(DisplRCP),each=dim(DisplRCPSummTbl[[1]])[1]),
                              do.call("rbind",DisplRCPSummTbl))
DisplRCPSummTbl$Year <- factor(DisplRCPSummTbl$Year,
                               levels=rev(seq(2099, 2299,by=50)))
## Divergence Novelty
DiveRCPSummTbl <- lapply(1: dim(DiveRCP)[3],
                         function(x){
                           out <- as.data.frame(table(factor(values(DiveRCP[[x]],na.rm=T),c(0,seq(2099, 2299,by=50))))/sum(table(values(DiveRCP[[x]],na.rm=T))))[-1,]
                           out$Freq <- round(out$Freq*100,2)
                           names(out) <- c("Year","Freq")
                           return(out)
                         })
DiveRCPSummTbl <- data.frame(RCP = rep(names(DiveRCP),each=dim(DiveRCPSummTbl[[1]])[1]),
                             do.call("rbind",DiveRCPSummTbl))
DiveRCPSummTbl$Year <- factor(DiveRCPSummTbl$Year,
                              levels=rev(seq(2099, 2299,by=50)))


#####
# Plot by when a given area will be novel
MDAllRCP2 <- resample(MDAllRCP,DisplRCP,method = "near")

AllRast <- c(MDAllRCP2,DisplRCP,DiveRCP)

names(AllRast) <- paste(rep(c("Climate Novelty","Dispacement","Dovergecne"),each=4),c("RCP26", "RCP45", "RCP60", "RCP85"),sep="_")
pdf("./Results/PDF/Sup_Figure_3.pdf",
    width = 10, height=10)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = as.factor(AllRast)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_manual(values = c("grey50",RColorBrewer::brewer.pal(5,"RdYlBu")),na.value = NA,#Do not map NA
                    breaks = c(0,seq(2099, 2299,by=50)), # Legend breaks
                    labels = c("No Change",seq(2099, 2299,by=50)), # Legend Labels
                    name = "Year when\nNo-Analogue?") + # Legend Title +
  facet_wrap(~lyr, ncol=3,nrow=4,dir="v") +
  geom_spatvector(fill = NA)
dev.off()