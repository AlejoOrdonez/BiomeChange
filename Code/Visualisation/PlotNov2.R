rm(list=ls());gc()
require(terra)
library(tidyverse)
library(tidyterra)
require(maptools)
data(wrld_simpl)
wrld_simpl2 <- vect(wrld_simpl)
wrld_simpl2 <- project(wrld_simpl2,"+proj=eck4")
setwd("/Volumes/MacPro 2013 Backup/BiomeChange/")
setwd("/Users/alejandroordonez/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")
Biome <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_100km.tif")

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



Model <- c("BIOCLIM","BIOCLIM2","koeppen_geiger","Seasonal")[4]#(Model <- c("BIOCLIM","BIOCLIM2","koeppen_geiger","Seasonal")[1])
####

# Estimate by when a given area will be novel
## Novelty estimated using a Mahalanobis distance 
MDAllRCP <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                   function(RCP){#(RCP <- "RCP60")
                     # Load the Treshold
                     Tresh <- readRDS(paste0("./Results2/Novelty/AllModels/",
                                             RCP,"/",Model,
                                             "/AllModels_",RCP,"_TreshSumm.rds"))
                     # Load the MD distances for each yeat                      
                     Tresh <- readRDS(paste0("./Results2/Novelty/AllModels/",
                                             RCP,"/",Model,
                                             "/AllModels_",RCP,"_TreshSumm.rds"))
                     
                     MDSum.ModlRCP <- lapply(paste0("./Results2/Novelty/AllModels/",
                                                    RCP,"/",Model,"/AllModels_",RCP,"_",
                                                    seq(2099, 2299,by=50),"_MDminSumm.tif"),
                                             function(x){rast(x,lyrs = "MD.min")})
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

# Plot by when a given area will be novel
pdf("/Volumes/MacPro 2013 Backup/BiomeChange/Results2/PDF/Fig1.pdf",
    width = 5, height=9)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = as.factor(MDAllRCP)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_d(palette = "muted", #color scheme 
                        na.value = NA,#Do not map NA
                        breaks = c(0,seq(2099, 2299,by=50)), # Legend breaks
                        labels = c("No Change",seq(2099, 2299,by=50)), # Legend Labels
                        name = "Year when\nNo-Analogue?" # Legend Title
  ) +
  facet_wrap(~lyr, ncol=1) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = Model) # Fig title
dev.off()
####
####
# Estimate how much of the land will be novel by 2300
colSums(apply(values(MDAllRCP,na.rm=T),
              2,
              function(x){out <- round(table(factor(x,c(0,seq(2099, 2299,by=50))))/sum(table(x))*100,2)})[-1,])

# Novel area by Yr and RCP as a table
MDAllRCPContTbl <- apply(values(MDAllRCP,na.rm=T),
                         2,
                         function(x){out <- round(table(factor(x,c(0,seq(2099, 2299,by=50))))/sum(table(x))*100,2)})[-1,]
# Rate of change
(MDAllRCPContTbl[-dim(MDAllRCPContTbl)[1],] - MDAllRCPContTbl[-1,])/50
mean(((MDAllRCPContTbl[-dim(MDAllRCPContTbl)[1],] - MDAllRCPContTbl[-1,])/50)[1:2,3:4])
mean(((MDAllRCPContTbl[-dim(MDAllRCPContTbl)[1],] - MDAllRCPContTbl[-1,])/50)[3:4,3:4])

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

# Avg yearly change for RCP 6.0 & 8.5 between 2010-2050
mean(tapply(MDAllRCPSummTbl$Freq/50,
       list(MDAllRCPSummTbl$RCP,MDAllRCPSummTbl$Year),
       mean)[1:2,4:5])

# Avg yearly change for RCP 6.0 & 8.5 between 2010-2050
mean(tapply(MDAllRCPSummTbl$Freq/50,
            list(MDAllRCPSummTbl$RCP,MDAllRCPSummTbl$Year),
            mean)[3:4,2:3])
####
####
# Estimate the percentage change of the area of each biome for each evaluate year
BiomePerChng <- lapply(MDAllRCP,
                       function(x){
                         BiomeSumm <- rbind(table(values(Biome,na.rm=T)), # Number of cells is a given biome
                                            (table(values(x,na.rm=T),values(Biome,na.rm=T)))[-1,]) # Changed cels per year in a given biome
                         # Estimate the prportion fo a given biome that changes
                         BiomeSumm <- apply(BiomeSumm,2,function(x){x[-1]/x[1]})
                         colnames(BiomeSumm) <- BiomeNames[,2]
                         return(BiomeSumm)
                         
                       })
names(BiomePerChng) <- c("RCP26", "RCP45", "RCP60", "RCP85")

# Summary per Biome of the total amount of change per 2300 for each Biome
BiomePerChngSumm <- sapply(BiomePerChng,
                           function(x){
                             round(apply(x,2,sum)*100,3)         
                           })
# Range of Biome change
apply(BiomePerChngSumm,2,range)



# Make a table of the biome proportion of novelty 
BiomePerChngTble <- lapply(BiomePerChng,
                           function(x){
                             data.frame(Biome = rep(1:14,each=dim(x)[1]),#rep(colnames(x),each=dim(x)[1]),
                                        Year = rep(rownames(x),dim(x)[2]),
                                        Freq = do.call("c",
                                                       lapply(1:dim(x)[2],
                                                              function(i){
                                                                x[,i]
                                                              })))
                           })
BiomePerChngTble <-  data.frame(RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"),
                                          each=dim(BiomePerChngTble[[1]])[1]),
                                do.call("rbind",BiomePerChngTble))
BiomePerChngTble$Year <- factor(BiomePerChngTble$Year,
                                rev(seq(2099, 2299,by=50)))
BiomePerChngTble$Freq <- round(BiomePerChngTble$Freq*100,2)

BiomePerChngTble <- rbind(data.frame(RCP = MDAllRCPSummTbl$RCP,
                                     Biome = rep(0,dim(MDAllRCPSummTbl)[1]),
                                     MDAllRCPSummTbl[,c("Year","Freq")]),
                          BiomePerChngTble)
BiomePerChngTble <- BiomePerChngTble[order(BiomePerChngTble$RCP),]
BiomePerChngTble$Biome <- BiomePerChngTble$Biome+1

pdf("/Volumes/MacPro 2013 Backup/BiomeChange/Results2/PDF/Fig2.pdf",
    width = 5, height=10)#width = 7, height=10)

ggplot(BiomePerChngTble, aes(fill=Year, y=Freq, x=Biome)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "RdYlBu") +
  ylab("Proportion Land considered novel") +
  xlab("Biome") +
  scale_x_discrete(limits= c("All",BiomeNames$Name)) +
  ggtitle(Model) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 2, unit = "cm")) +
  facet_wrap(~RCP,ncol=1)
dev.off()



a <- (BiomePerChngTble[BiomePerChngTble$RCP=="RCP26",c("Year","Perc")])
a$Year <- as.numeric(as.character(a$Year))
plot(a)

abline(lm(a[,2:1]))
(mean(a$Perc/50))


range(apply(tapply(MDAllRCPSummTbl$Freq/50,
             list(MDAllRCPSummTbl$RCP,MDAllRCPSummTbl$Year),
             mean),2,mean))

# Avg yearly change for RCP 2.6 & 4.5 between 2010-2050
mean(tapply(MDAllRCPSummTbl$Freq/50,
            list(MDAllRCPSummTbl$RCP,MDAllRCPSummTbl$Year),
            mean)[1:2,4:5])
# Avg yearly change for RCP 6.0 & 8.5 between 2010-2050
mean(tapply(MDAllRCPSummTbl$Freq/50,
            list(MDAllRCPSummTbl$RCP,MDAllRCPSummTbl$Year),
            mean)[1:2,4:5])

# Avg yearly change for RCP 6.0 & 8.5 between 2010-2050
mean(tapply(MDAllRCPSummTbl$Freq/50,
            list(MDAllRCPSummTbl$RCP,MDAllRCPSummTbl$Year),
            mean)[3:4,2:3])

####




