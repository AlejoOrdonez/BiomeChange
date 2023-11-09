rm(list=ls());gc()
require(terra)
library(tidyterra)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")

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
                       #Load the divergecne for each year
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

# Estimate the percentage change of the area of each biome for each evaluate year
# Load the biome data
# Load the Biome rast and shape
Biome <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_100km.tif")
Biome <- resample(Biome,MDAllRCP,method="near")
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
## Estimate biome proportions for Climate novelty

BiomePerChng <- lapply(1:4,
                       function(x){
                         BiomeSumm <- rbind(table(values(Biome)), # Number of cells is a given biome
                                            table(factor(values(MDAllRCP[[x]]),levels=c(0,seq(2099, 2299,by=50))),values(Biome))[-1,]) # Changed cels per year in a given biome
                         # Estimate the prportion fo a given biome that changes
                         BiomeSumm <- apply(BiomeSumm,2,function(x){x[-1]/x[1]})
                         colnames(BiomeSumm) <- BiomeNames[,2]
                         BiomeTbl <- data.frame(Biome = rep(1:14,each=dim(BiomeSumm)[1]),
                                                Year = rep(rownames(BiomeSumm),dim(BiomeSumm)[2]),
                                                Freq = do.call("c",
                                                               lapply(1:dim(BiomeSumm)[2],
                                                                      function(i){
                                                                        BiomeSumm[,i]
                                                                      })))
                         return(BiomeTbl)
                       })
names(BiomePerChng) <- c("RCP26", "RCP45", "RCP60", "RCP85")
BiomePerChngTble <-  data.frame(RCP = sapply(strsplit(row.names(do.call("rbind",BiomePerChng)),"[.]"),function(x){x[1]}),
                                do.call("rbind",BiomePerChng))
BiomePerChngTble$Year <- factor(BiomePerChngTble$Year,
                                rev(seq(2099, 2299,by=50)))
BiomePerChngTble$Freq <- round(BiomePerChngTble$Freq*100,2)

# Merge the All sites and Biome table
MDPerChngTble <- rbind(data.frame(RCP = MDAllRCPSummTbl$RCP,
                                  Biome = rep(0,dim(MDAllRCPSummTbl)[1]),
                                  MDAllRCPSummTbl[,c("Year","Freq")]),
                       BiomePerChngTble)
MDPerChngTble <- MDPerChngTble[order(MDPerChngTble$RCP),]
MDPerChngTble$Biome <- MDPerChngTble$Biome+1
MDPerChngTble$Method <- "Climate Novelty"

## Estimate biome proportions for Displacement 
Biome <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_100km.tif")
Biome <- resample(Biome,DisplRCP,method="near")

BiomePerChng <- lapply(1:4,
                       function(x){
                         
                         BiomeSumm <- rbind(table(values(Biome)), # Number of cells is a given biome
                                            table(factor(values(DisplRCP[[x]]),levels=c(0,seq(2099, 2299,by=50))),values(Biome))[-1,]) # Changed cels per year in a given biome
                         # Estimate the prportion fo a given biome that changes
                         BiomeSumm <- apply(BiomeSumm,2,function(x){x[-1]/x[1]})
                         colnames(BiomeSumm) <- BiomeNames[,2]
                         BiomeTbl <- data.frame(Biome = rep(1:14,each=dim(BiomeSumm)[1]),
                                                Year = rep(rownames(BiomeSumm),dim(BiomeSumm)[2]),
                                                Freq = do.call("c",
                                                               lapply(1:dim(BiomeSumm)[2],
                                                                      function(i){
                                                                        BiomeSumm[,i]
                                                                      })))
                         return(BiomeTbl)
                       })
names(BiomePerChng) <- c("RCP26", "RCP45", "RCP60", "RCP85")
BiomePerChngTble <-  data.frame(RCP = sapply(strsplit(row.names(do.call("rbind",BiomePerChng)),"[.]"),function(x){x[1]}),
                                do.call("rbind",BiomePerChng))
BiomePerChngTble$Year <- factor(BiomePerChngTble$Year,
                                rev(seq(2099, 2299,by=50)))
BiomePerChngTble$Freq <- round(BiomePerChngTble$Freq*100,2)

# Merge the All sites and Biome table
DisplPerChngTble <- rbind(data.frame(RCP = DisplRCPSummTbl$RCP,
                                     Biome = rep(0,dim(DisplRCPSummTbl)[1]),
                                     DisplRCPSummTbl[,c("Year","Freq")]),
                          BiomePerChngTble)
DisplPerChngTble <- DisplPerChngTble[order(DisplPerChngTble$RCP),]
DisplPerChngTble$Biome <- DisplPerChngTble$Biome+1
DisplPerChngTble$Method <- "Displacement"

## Estimate biome proportions for Divergence 
Biome <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_100km.tif")
Biome <- resample(Biome,DiveRCP,method="near")

BiomePerChng <- lapply(1:4,
                       function(x){
                         
                         BiomeSumm <- rbind(table(values(Biome)), # Number of cells is a given biome
                                            table(factor(values(DiveRCP[[x]]),levels=c(0,seq(2099, 2299,by=50))),values(Biome))[-1,]) # Changed cels per year in a given biome
                         # Estimate the prportion fo a given biome that changes
                         BiomeSumm <- apply(BiomeSumm,2,function(x){x[-1]/x[1]})
                         colnames(BiomeSumm) <- BiomeNames[,2]
                         BiomeTbl <- data.frame(Biome = rep(1:14,each=dim(BiomeSumm)[1]),
                                                Year = rep(rownames(BiomeSumm),dim(BiomeSumm)[2]),
                                                Freq = do.call("c",
                                                               lapply(1:dim(BiomeSumm)[2],
                                                                      function(i){
                                                                        BiomeSumm[,i]
                                                                      })))
                         return(BiomeTbl)
                       })
names(BiomePerChng) <- c("RCP26", "RCP45", "RCP60", "RCP85")
BiomePerChngTble <-  data.frame(RCP = sapply(strsplit(row.names(do.call("rbind",BiomePerChng)),"[.]"),function(x){x[1]}),
                                do.call("rbind",BiomePerChng))
BiomePerChngTble$Year <- factor(BiomePerChngTble$Year,
                                rev(seq(2099, 2299,by=50)))
BiomePerChngTble$Freq <- round(BiomePerChngTble$Freq*100,2)

# Merge the All sites and Biome table
DivePerChngTble <- rbind(data.frame(RCP = DiveRCPSummTbl$RCP,
                                    Biome = rep(0,dim(DiveRCPSummTbl)[1]),
                                    DiveRCPSummTbl[,c("Year","Freq")]),
                         BiomePerChngTble)
DivePerChngTble <- DivePerChngTble[order(DivePerChngTble$RCP),]
DivePerChngTble$Biome <- DivePerChngTble$Biome+1
DivePerChngTble$Method <- "Divergence"

## PLot the ampount per year
AlPerChngTblel <- rbind(MDPerChngTble,
                        DisplPerChngTble,
                        DivePerChngTble)

pdf("./Results/PDF/Sup_Figure_2.pdf",
    width = 10, height=12)#width = 10, height=5)
ggplot(AlPerChngTblel, aes(fill=Year, y = Freq, x=Biome)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~RCP+Method, ncol=3) +
  scale_x_discrete(limits= c("All",BiomeNames$Name)) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) + 
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(5,"RdYlBu")),na.value = NA,#Do not map NA
                    breaks = c(seq(2099, 2299,by=50)), # Legend breaks
                    labels = c(seq(2099, 2299,by=50)), # Legend Labels
                    name = "Year when\nNo-Analogue?") + # Legend Title
  ylab("Proportion Land considered novel") +
  xlab("Biome") 
dev.off()