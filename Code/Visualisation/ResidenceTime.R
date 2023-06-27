rm(list=ls());gc()
require(terra)
library(tidyterra)
library(tidyverse)

# Load a simpified map of the world
require(maptools)
# Load the Biome rast and shape
setwd("/Volumes/MacPro 2013 Backup/BiomeChange")
BiomeShp <-vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")
BiomeShp <- BiomeShp[c(BiomeShp$BIOME<15),]
BiomeShp <- project(BiomeShp,"+proj=eck4")

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

## Displacment 
Model <- c("BIOCLIM","BIOCLIM2","koeppen_geiger","Seasonal")[4]#(Model <- c("BIOCLIM","BIOCLIM2","koeppen_geiger","Seasonal")[1])
DisplARM <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                   function(RCP){#(RCP <- "RCP26")
                     #Load the displacement for each year
                     RastTmpList <- lapply(paste0("./Results2/Velocity/AllModels/",RCP,"/",Model,
                                                  "/AllModels_DispDiv_",RCP,"_",c(seq(2049, 2299,by=50)-50),"to",seq(2049, 2299,by=50),".tif"),
                                           function(x){rast(x, lyrs="Displ.ARM")})
                     RastTmp <- do.call("c",RastTmpList)
                     names(RastTmp) <- paste(c(seq(2049, 2299,by=50)-50),"to",seq(2049, 2299,by=50))
                     return(RastTmp)})
names(DisplARM) <- c("RCP26", "RCP45", "RCP60", "RCP85")

if(!paste0("ResTimebyRCP_",Model,".rds") %in%dir("./Results2/Velocity/AllModels/")){
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
          paste0("./Results2/Velocity/AllModels/ResTimebyRCP_",Model,".rds"))
} else{
  ResTimebyRCP <- readRDS(paste0("./Results2/Velocity/AllModels/ResTimebyRCP_",Model,".rds"))
}

# Box plot of mean displacement per Ecoregion
BiomeShp2 <- BiomeShp[ResTimebyRCP[[1]]$ID,] # prun the safile to only polygons with displacment data

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
  ggplot(ResTimebyRCPTbl, aes(x=Mean, y=Biome)) + 
  geom_boxplot() +
  scale_x_continuous(trans='log10',
                     name = "Residence time in Years",
                     breaks = c(0,1,10,100,1000)) +
    facet_wrap("RCP")

  
  # Proportion of ecoregions with residence time in the order of Decades
  sapply(ResTimebyRCP,
         function(x){
           SummTmp <- tapply(x$Mean,BiomeShp2$ECO_ID[x$ID],mean)
           round(100*sum(round(log10(SummTmp))==1)/length(SummTmp),2)
         })
  #  Proportion of ecoregions with residence time above centenial
  sapply(ResTimebyRCP,
         function(x){
           SummTmp <- tapply(x$Mean,BiomeShp2$ECO_ID[x$ID],mean)
           round(100*sum(round(log10(SummTmp))>1)/length(SummTmp),2)
         })
  sapply(ResTimebyRCP,
         function(x){#(x<-ResTimebyRCP[[4]])
           sum(round(100*table(round(log10(x$Mean)))/dim(x)[1],2)[as.character(2:4)],na.rm=T)
         })
  # Proportion of ecoregions with residence time in given time category
  lapply(ResTimebyRCP,
         function(x){
           SummTmp <- tapply(x$Mean,BiomeShp2$ECO_ID[x$ID],mean)
           round(100*(table(round(log10(SummTmp)))/length(SummTmp)),2)
         })
           
  # Proportion of ecoregions ina given biome with residence time in given time category
  lapply(ResTimebyRCP,
         function(x){
           SummTmp <- tapply(x$Mean,BiomeShp2$ECO_ID[x$ID],mean)
           out <- table(round(log10(SummTmp)),BiomeShp2$BIOME[(match(as.numeric(names(SummTmp)),BiomeShp2$ECO_ID))])
           round(100*(apply(out,2,function(x){x/sum(x)})),2)
         })
  
  # pdf("~/Desktop/Fig.pdf")
  # a <- round(log10(tapply(ResTimebyRCP[[4]]$Mean,BiomeShp2$ECO_ID,mean)))
  # BiomeShp2$Mean <-factor(a,
  #                         levels=-3:3)
  # 
  # ggplot(BiomeShp2) +
  #   geom_spatvector(aes(fill = Mean),col=NA) +
  #   scale_fill_viridis_d(breaks=-3:3)
  # dev.off()  

  