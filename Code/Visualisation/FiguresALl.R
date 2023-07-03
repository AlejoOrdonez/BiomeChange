rm(list=ls());gc()
require(terra)
library(tidyverse)
library(tidyterra)
require(maptools)
data(wrld_simpl)
wrld_simpl2 <- vect(wrld_simpl)
wrld_simpl2 <- project(wrld_simpl2,"+proj=eck4")
#setwd("/Volumes/MacPro 2013 Backup/BiomeChange/")
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange")

# Load the Biome rast and shape
Biome <- rast("./Data/WWF-Biomes/WWF_BIOME_eck4_100km.tif")
BiomeShp <-vect("./Data/WWF-Biomes/wwf_terr_ecos.shp")
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

#####
#Figure 1
# Plot by when a given area will be novel
pdf("./Results2/PDF/Fig1.pdf",width = 5, height=10)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = as.factor(MDAllRCP)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_manual(values = c("grey50",RColorBrewer::brewer.pal(5,"RdYlBu")),na.value = NA,#Do not map NA
                    breaks = c(0,seq(2099, 2299,by=50)), # Legend breaks
                    labels = c("No Change",seq(2099, 2299,by=50)), # Legend Labels
                    name = "Year when\nNo-Analogue?") + # Legend Title +
  facet_wrap(~lyr, ncol=1) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = Model) # Fig title
dev.off()
#####
##### Figure 2
# Estimate the percentage change of the area of all areas for each evaluate year
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

# Estimate the percentage change of the area of each biome for each evaluate year
BiomePerChng <- lapply(MDAllRCP,
                       function(x){
                         BiomeSumm <- rbind(table(values(Biome,na.rm=T)), # Number of cells is a given biome
                                            (table(values(x,na.rm=T),values(Biome,na.rm=T)))[-1,]) # Changed cels per year in a given biome
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
BiomePerChngTble <-  data.frame(RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"),
                                          each=dim(BiomePerChng[[1]])[1]),
                                do.call("rbind",BiomePerChng))
BiomePerChngTble$Year <- factor(BiomePerChngTble$Year,
                                rev(seq(2099, 2299,by=50)))
BiomePerChngTble$Freq <- round(BiomePerChngTble$Freq*100,2)

# Merge the All sites and Biome table
PerChngTble <- rbind(data.frame(RCP = MDAllRCPSummTbl$RCP,
                                Biome = rep(0,dim(MDAllRCPSummTbl)[1]),
                                MDAllRCPSummTbl[,c("Year","Freq")]),
                          BiomePerChngTble)
PerChngTble <- PerChngTble[order(PerChngTble$RCP),]
PerChngTble$Biome <- PerChngTble$Biome+1

pdf("./Results2/PDF/Fig2.pdf",width = 5, height=13)#width = 7, height=10)
ggplot(PerChngTble, aes(fill=Year, y = Freq, x=Biome)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~RCP,ncol=1) +
  scale_x_discrete(limits= c("All",BiomeNames$Name)) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) + 
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(5,"RdYlBu")),na.value = NA,#Do not map NA
                    breaks = c(seq(2099, 2299,by=50)), # Legend breaks
                    labels = c(seq(2099, 2299,by=50)), # Legend Labels
                    name = "Year when\nNo-Analogue?") + # Legend Title
  ylab("Proportion Land considered novel") +
  xlab("Biome") +
  ggtitle(Model)
dev.off()

#####
#####
rm(list= ls()[!ls()%in%c("wrld_simpl2","Biome","BiomeShp","BiomeNames","Model")]);gc()
# Plot the mean displacemnt estimated over the evalyated 50yrs periods 
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

DisplARMList <- lapply(DisplARM,function(x){mean(x)})
DisplARMMeanRast <- c(DisplARMList[[1]],
                      DisplARMList[[2]],
                      DisplARMList[[3]],
                      DisplARMList[[4]])
names(DisplARMMeanRast) <- c("RCP26", "RCP45", "RCP60", "RCP85")

pdf("./Results2/PDF/Fig3.pdf",width = 5, height=10)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = log10(DisplARMMeanRast)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = -2:2, # Legend breaks
                        labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Speed [km/yrs]", # Legend Title
                        limit = c(-2,2)) +
  facet_wrap(~lyr,ncol=1) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  theme(legend.position="bottom") +
  labs(title = Model) # Fig title
dev.off()
#####
#####
rm(list= ls()[!ls()%in%c("wrld_simpl2","Biome","BiomeShp","BiomeNames","Model","DisplARM")]);gc()
# Plot the Mean velocity per time period
## Summary of Displacements globaly
DisARMSumm <- lapply(DisplARM,
                     function(x){
                       out<- data.frame(Year = factor(names(x)),
                                        DisplArm = 10^apply(log10(values(x,na.rm=T)),
                                                            2,
                                                            mean))
                     })
### Compile geometric means of the displacement per Yr/RCP as a dataframe
DisARMSumm <- do.call("rbind",DisARMSumm)
DisARMSumm<- data.frame(RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"),
                                   each=dim(DisARMSumm)[1]/4),
                        DisARMSumm)

## Summary of Displacements per biome
DisARMBiomeSumm <- lapply(DisplARM,
                          function(x){
                            out <- lapply(1:6,
                                          function(i){#(i<-1)
                                            10^tapply(log10(values(x[[i]])),
                                                      values(Biome),
                                                      mean,na.rm=T)[-15]
                                          })
                            out <- data.frame(Biome = BiomeNames$Name,
                                              Year = rep(names(x),each=14),
                                              DisplArm = do.call("c",out))
                            return(out)})

# merge Global and per biome assessments
DisARMAll <- rbind(data.frame(RCP = DisARMSumm$RCP,
                              Biome = "All",
                              DisARMSumm[,-1]),
                   data.frame (RCP = rep(names(DisplARM), each = dim(DisARMBiomeSumm[[1]])[1]),
                               do.call("rbind",DisARMBiomeSumm)))
DisARMAll$Biome <- factor(DisARMAll$Biome)
DisARMAll <- DisARMAll[order(DisARMAll$Biome,decreasing = T),]


# plot the Rate of Displacement 
pdf("./Results2/PDF/Fig4.pdf",
    width = 7.5, height=10)#width = 10.5, height=7)
ggplot(DisARMAll,
       aes(x=Year, y=DisplArm, group=Biome ,color = Biome)) + 
  ylab("Displaecment\n[km/year]") +
  scale_color_manual(values=c("Black",scales::brewer_pal(palette = "Dark2")(8),scales::brewer_pal(palette = "Paired")(6))) +
  
  scales::brewer_pal(palette = "Dark2")
  xlab("Time Period") +
  geom_line() +
  geom_point(shape = 19) +
  scale_y_continuous(trans='log10') + 
  #scale_color_manual(values=c("darkgreen","purple","gold","red")) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_wrap(facets = "RCP",ncol=1) +
  labs(title = Model) # Fig title
dev.off()
#####
#####
rm(list= ls()[!ls()%in%c("wrld_simpl2","Biome","BiomeShp","BiomeNames","Model","DisplARM")]);gc()
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
BiomeShp2 <- BiomeShp[ResTimebyRCP[[1]]$ID,] # prune the shapefile to only polygons with displacement data
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
pdf("./Results2/PDF/Fig5.pdf",
    width = 9, height=8)
ggplot(ResTimebyRCPTbl, aes(x=Mean, y=Biome)) + 
  geom_boxplot() +
  scale_x_continuous(trans='log10',
                     name = "Residence time in Years",
                     breaks = c(0,1,10,100,1000)) +
  facet_wrap("RCP") +
  labs(title = Model) # Fig title
dev.off()
#####
#####
rm(list= ls()[!ls()%in%c("wrld_simpl2","Biome","BiomeShp","BiomeNames","Model")]);gc()
# Maps of divergence angles
Divergence <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                     function(RCP){#(RCP <- "RCP26")
                       #Load the displacement for each year
                       RastTmpList <- lapply(paste0("./Results2/Velocity/AllModels/",RCP,"/",Model,
                                                    "/AllModels_DispDiv_",RCP,"_",c(seq(2049, 2299,by=50)-50),"to",seq(2049, 2299,by=50),".tif"),
                                             function(x){rast(x, lyrs="Dive")})
                       RastTmp <- do.call("c",RastTmpList)
                       names(RastTmp) <- paste(c(seq(2049, 2299,by=50)-50),"to",seq(2049, 2299,by=50))
                       return(RastTmp)})
names(Divergence) <- c("RCP26", "RCP45", "RCP60", "RCP85")

# Plot the divergence
# Map the projected displacement
DivergenceSumm <- c(mean(Divergence[[1]]),
                    mean(Divergence[[2]]),
                    mean(Divergence[[3]]),
                    mean(Divergence[[4]]))
names(DivergenceSumm) <- names(Divergence)

# Map diveregences
pdf("./Results2/PDF/Fig6.pdf",
    width = 5, height=10)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = DivergenceSumm) + # Map the Displacement
  facet_wrap(~lyr,ncol = 1,dir= "v") +
  scale_fill_gradientn(colours = colorRampPalette(c("blue","red","yellow"))(255),
                      na.value = NA,
                        breaks = c(0,60,120,180), # Legend breaks
                        #labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Divergence angle", # Legend Title
                        limit = c(0,180)) +
  theme(legend.position="bottom") +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = Model) # Fig title
dev.off()
#####
#####
rm(list= ls()[!ls()%in%c("wrld_simpl2","Biome","BiomeShp","BiomeNames","Model","DivergenceSumm")]);gc()
# Proportion of cells with Congruent (0-60), orthogonal (60-120), and opposing (120-180) divergences per biome
DivSumByBiome <- lapply(DivergenceSumm,
                        function(x){#(x <- DivergenceSumm[[1]])
                          a <- terra::extract(x, BiomeShp, fun = mean,na.rm=T)
                          a <- tapply(a[,2],BiomeShp$ECO_ID,mean,na.rm=T)
                          # Summary All sites
                          y <- factor(ifelse(a<60,1,ifelse(a<120,2,3)),
                                      levels = 1:3)
                          y <- as.numeric(round(100*(table(y)/length(y)),2))
                          out <- data.frame(Biome = "All",
                                            Sind = c("Congruent","Orthogonal","Opposing"),
                                            Freq = y) 
                          # Summary by Biome
                          y2 <- lapply(1:14,
                                       function(i){#(i<-2)
                                         j <- a[BiomeShp$BIOME[match(names(a),BiomeShp$ECO_ID)]==i]
                                         j <- factor(ifelse(j<60,1,ifelse(j<120,2,3)),
                                                     levels = 1:3)
                                         j <- as.numeric(round(100*(table(j)/length(j)),2))
                                         j <- data.frame(Biome = BiomeNames$Name[i],
                                                         Sind = c("Congruent","Orthogonal","Opposing"),
                                                         Freq = j) 
                                         return(j)
                                       })
                          out <- rbind(out,do.call("rbind",y2))
                          return(out)})
DivSumByBiome2 <- data.frame(RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"),
                                       each = dim(DivSumByBiome[[1]])[1]),
                             do.call("rbind",DivSumByBiome))

pdf("./Results2/PDF/Fig7.pdf",
    width = 5, height=13)
ggplot(DivSumByBiome2, aes(fill=Sind, y=Freq, x=Biome)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~RCP,ncol=1) +
  scale_fill_manual(values = c("blue","yellow","red"),
                    na.value = NA,#Do not map NA
                    #labels = c("Congruent","Orthogonal","Opposing"), # Legend Labels
                    name = "Divergecne class") + # Legend Title +
  ylab("Proportion Land under a divergecne syndrome") +
  xlab("Time period") +

  theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  ggtitle(Model)
dev.off()
#####
#####
# Estimate by when a given area will be novel - Novelty estimated using a Mahalanobis distance 
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
# Turn into novel area (Novel at any point)
MDAllRCP <- (MDAllRCP!=0)*1


# Estimate the mean Displacement over the evaluated time periods (2049,2099,2149,2199,2249,2299) for all RCPs - Displacement estimated using an ARM model
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
                       RastTmpList <- lapply(paste0("./Results2/Velocity/AllModels/",RCP,"/",Model,
                                                    "/AllModels_DispDiv_",RCP,"_",c(seq(2049, 2299,by=50)-50),"to",seq(2049, 2299,by=50),".tif"),
                                             function(x){rast(x, lyrs="Dive")})
                       RastTmp <- do.call("c",RastTmpList)
                       names(RastTmp) <- paste(c(seq(2049, 2299,by=50)-50),"to",seq(2049, 2299,by=50))
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

NovelCriteriaList <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                            function(x){
                              c(MDAllRCP[[x]],
                                DisplARMMeanRast[[x]],
                                DivergenceSumm[[x]])
                            })
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

pdf("./Results2/PDF/Fig8.pdf",
    width = 5, height=13)#width = 10, height=5)
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
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = Model) # Fig title
dev.off()
#####
#####
# Proportion of cell per mechanism _ALL
=
#####
