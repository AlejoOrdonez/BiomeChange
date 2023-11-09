rm(list=ls());gc()
require(terra)
library(tidyterra)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")

# Load the Biome rast and shape
#Biome <- rast("WWF-Biomes/WWF_BIOME_eck4_100km.tif")
BiomeShp <- vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")
BiomeShp <- project(BiomeShp,"+proj=eck4")
BiomeShp <- BiomeShp[BiomeShp$BIOME<15]
BiomeNames <- data.frame(ID = 1:14,
                         Name = c("Tropical & Subtropical Moist Broadleaf Forests",
                                  "Tropical & Subtropical Dry Broadleaf Forests",
                                  "Tropical & Subtropical Coniferous Forests",
                                  "Temperate Broadleaf & Mixed Forests",
                                  "Temperate Conifer Forests",
                                  "Boreal Forests/Taiga",
                                  "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                  "Temperate Grasslands, Savannas & Shrublands",
                                  "Flooded Grasslands & Savannas",
                                  "Montane Grasslands & Shrublands",
                                  "Tundra",
                                  "Mediterranean Forests, Woodlands & Scrub",
                                  "Deserts & Xeric Shrublands",
                                  "Mangroves"))

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
names(DisplARM) <- c("RCP26", "RCP45", "RCP60", "RCP85")

if(!paste0("ResTimebyRCP_Seasonal.rds") %in%dir("./Results/Displacement_Divergence/Mean_All_Models/")){
  ResTimebyRCP <- lapply(DisplARM,
                         function(x){#(x<-DisplARM[[4]])
                           ResTimeList <- lapply(1:dim(BiomeShp)[1],function(i){#(i<-c(1:dim(BiomeShp)[1])[1])
                                                   BiomeTmp <- BiomeShp[i,]
                                                   DispPerBiome <- apply(terra::extract(x,BiomeTmp)[,-1],2,median,na.rm=T)
                                                   #Radius <- mean(sqrt(BiomeTmp$area_km2/pi))
                                                   Radius <- mean(dist(crds(BiomeTmp)))/1000
                                                   ResTime <- Radius/DispPerBiome
                                                   return(ResTime)
                                                 })
                           ResTimeTbl <- data.frame(ID =1:dim(BiomeShp)[1],
                                                    Biome = factor(BiomeNames$Name[BiomeShp$BIOME]),
                                                    do.call("rbind",ResTimeList))
                           ResTimeTbl$Mean <- apply(ResTimeTbl[,-c(1:2)],1,mean,na.rm=T)
                           ResTimeTbl <- na.omit(ResTimeTbl)
                           return(ResTimeTbl)
                         })
  saveRDS(ResTimebyRCP,
          paste0("./Results/Displacement_Divergence/Mean_All_Models/ResTimebyRCP_Seasonal.rds"))
} else{
  ResTimebyRCP <- readRDS(paste0("./Results/Displacement_Divergence/Mean_All_Models/ResTimebyRCP_Seasonal.rds"))
}

# Box plot of mean displacement per Ecoregion
BiomeShp2 <- BiomeShp[ResTimebyRCP[[1]]$ID,] # prune the saved file to only polygons with displacement data

#Summary of displacement by Ecoregion/Biome
ResTimebyRCPList <- lapply(ResTimebyRCP,
                           function(x){
                             ResTimebyRCPTmp <- tapply(x$Mean,BiomeShp2$ECO_ID,mean)
                             ResTimebyRCPTmp <- data.frame(Biome = factor(BiomeNames$Name[BiomeShp2$BIOME[match(names(ResTimebyRCPTmp),BiomeShp2$ECO_ID)]]),
                                                           Mean = as.numeric(ResTimebyRCPTmp))
                             return(ResTimebyRCPTmp)                             
                           })

ResTimebyRCPTbl <- data.frame(RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"),
                                        each = dim(ResTimebyRCPList[[1]])[1]),
                              do.call("rbind",ResTimebyRCPList))

# plot the distribution of displacements
pdf("./Results/PDF/Figure_3.pdf",
    width = 9, height=8)
ggplot(ResTimebyRCPTbl, aes(x=Mean, y=Biome)) + 
  geom_boxplot() +
  scale_x_continuous(trans='log10',
                     name = "Residence time in Years",
                     breaks = c(0,1,10,100,1000)) +
  facet_wrap("RCP")
dev.off()
  