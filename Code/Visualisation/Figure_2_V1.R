rm(list=ls());gc()
require(terra)
library(tidyterra)
library(tidyverse)

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
                     # cumulative No nalogue
                     NoAnalogue_comulative <- app(NoAnalogue,
                                                  function(x){
                                                    out <- x
                                                    if(!is.na(x[1]) & sum(x)>0){
                                                      out[min(which(x==1)):5] <- 1
                                                    }
                                                    return(out)
                                                  })
                     # Proportion of no analoge for a given year
                     NoAnalYr <- apply(values(NoAnalogue_comulative,na.rm=T),2,sum)/dim(values(NoAnalogue,na.rm=T))[1]
                     return(NoAnalYr)
                   })
# <Mean per year rates of change - avg of change per 50yr bins 
lapply(MDAllRCP,function(x){mean((x[-1]-x[-5])*100)/50})
# Make the summary
MDAllRCP <- data.frame(Mechanism = "Climate Novelty",
                       RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"), each=5),
                       Year = as.numeric(rep(names(MDAllRCP[[1]]), 4)),
                       Proportion = as.numeric(unlist(MDAllRCP))
                       )


ggplot(MDAllRCP) +aes(x=Year,y=Proportion,colour =RCP)+geom_smooth()
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
                     names(RastTmp) <- seq(2099, 2299,by=50)+1
                     # Assess if the Velocity corssess the treshold
                     NoAnalogue <- RastTmp > 1
                     # cumulative No nalogue
                     NoAnalogue_comulative <- app(NoAnalogue,
                                                  function(x){
                                                    out <- x
                                                    if(!is.na(x[1]) & sum(x)>0){
                                                      out[min(which(x==1)):5] <- 1
                                                    }
                                                    return(out)
                                                  })
                     # Proportion of no analoge for a given year
                     NoAnalYr <- apply(values(NoAnalogue_comulative,na.rm=T),2,sum)/dim(values(NoAnalogue,na.rm=T))[1]
                      return(NoAnalYr)
                     })
names(DisplARM) <- c("RCP26", "RCP45", "RCP60", "RCP85")
# <Mean per year rates of change - avg of change per 50yr bins 
lapply(DisplARM,function(x){mean((x[-1]-x[-5])*100)/50})

# Make the summary
DisplAllRCP <- data.frame(Mechanism = "Displacement",
                       RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"), each=5),
                       Year = as.numeric(rep(names(DisplARM[[1]]), 4)),
                       Proportion = as.numeric(unlist(DisplARM)))
ggplot(DisplAllRCP) +aes(x=Year,y=Proportion,colour =RCP)+geom_smooth()
# Estimate the mean divergence accorss evaluated time periods (2049,2099,2149,2199,2249,2299) for all RCPs
# Displacement estimated using an ARM model

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
                       # Assess if the Velocity corssess the treshold
                       NoAnalogue <- (RastTmp > 60 & RastTmp < 120)*1
                       # cumulative No nalogue
                       NoAnalogue_comulative <- app(NoAnalogue,
                                                    function(x){
                                                      out <- x
                                                      if(!is.na(x[1]) & sum(x)>0){
                                                        out[min(which(x==1)):5] <- 1
                                                      }
                                                      return(out)
                                                    })
                       # Proportion of no analoge for a given year
                       NoAnalYr <- apply(values(NoAnalogue_comulative,na.rm=T),2,sum)/dim(values(NoAnalogue,na.rm=T))[1]
                       return(NoAnalYr)
                       })
names(Divergence) <- c("RCP26", "RCP45", "RCP60", "RCP85")
# <Mean per year rates of change - avg of change per 50yr bins 
lapply(Divergence,function(x){mean((x[-1]-x[-5])*100)/50})
# Make the summary
DivAllRCP <- data.frame(Mechanism = "Divergence",
                          RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"), each=5),
                          Year = as.numeric(rep(names(Divergence[[1]]), 4)),
                          Proportion = as.numeric(unlist(Divergence)))

ggplot(DivAllRCP) +aes(x=Year,y=Proportion,colour =RCP)+geom_smooth()


pdf("./Results/PDF/Figure_2.pdf",
    width = 5, height=7)#width = 10, height=5)

AlltrendsRCP <- rbind(MDAllRCP,DisplAllRCP,DivAllRCP)

AlltrendsRCP$Proportion <- round(AlltrendsRCP$Proportion *100,2)
ggplot(AlltrendsRCP) +
  aes(x=Year,y=Proportion,col =Mechanism) +
  scale_colour_manual(values=c("#e41a1c", "#377eb8","#4daf4a")) +
  geom_smooth() +
  facet_wrap(facets= vars(RCP),ncol =1)

dev.off()
