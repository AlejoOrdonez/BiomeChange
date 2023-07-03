rm(list=ls());gc()
require(terra)
library(tidyterra)
library(tidyverse)
require(maptools)
data(wrld_simpl)
wrld_simpl2 <- vect(wrld_simpl)
wrld_simpl2 <- project(wrld_simpl2,"+proj=eck4")
rm(wrld_simpl);gc()
# Load the Biome rast and shape
setwd("/Volumes/MacPro 2013 Backup/BiomeChange")
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

pdf("/Volumes/MacPro 2013 Backup/BiomeChange/Results2/PDF/Fig8.pdf",
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
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = Model) # Fig title
dev.off()



# Proportion of cell per mechanism _ALL
NovSullAll <- apply(values(PlotTest),2,
                              function(x){round(100*(table(factor(x,levels = 1:8))/sum(table(x))),2)})

# Proportion of cell per mechanism Biome

NovSullBiome <- lapply(PlotTestList,
                       function(x){#(x<- PlotTestList[[4]])
                         out <- terra::extract(x, BiomeShp, fun = max,na.rm=T)
                         out2 <- apply(table(factor(out$lyr.1,levels=1:8),
                                              BiomeShp$BIOME),2,
                                       function(x){
                                         round(100*(x/sum(x)),2)
                                       })
                         data.frame(Biome = rep(1:14,each=dim(out2)[1]),
                                    Mech = 1:8,
                                    Frec = c(out2))
                       })

NovSull <- rbind(data.frame(RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"), each=8),
                            Biome = "All",
                            Mech = 1:8,
                            Frec = c(NovSullAll)),
                 data.frame(RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"), each=dim(NovSullBiome[[1]])[1]),
                            do.call("rbind",NovSullBiome))
                 )

NovSull$Biome[NovSull$Biome!="All"] <- BiomeNames$Name[as.numeric(NovSull$Biome[NovSull$Biome!="All"])]
NovSull$Biome <- factor(NovSull$Biome)
NovSull$Mech <- factor(NovSull$Mech,
                       levels = 1:8,
                       labels = c("No Novelty",
                                  "NovClim",
                                  "FastDis",
                                  "OrthDiv",
                                  "NovClim/FastDis",
                                  "FastDis/OrthDiv",
                                  "NovClim/OrthDiv",
                                  "NovClim/FastDis/OrthDiv"))

pdf("/Volumes/MacPro 2013 Backup/BiomeChange/Results2/PDF/Fig9.pdf",
    width = 5, height=10)#width = 7, height=10)

ggplot(NovSull, aes(fill=Mech, y=Frec, x=Biome)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_brewer(palette = "RdYlBu") +
  scale_fill_manual(values = c("#31a354","#54278f","#756bb1","#cbc9e2","#ffffb2","#ffeda0","#feb24c","#e41a1c"))+
  ylab("Proportion Land novel under a given mechanism") +
  xlab("Biome") +
  ggtitle(Model) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  facet_wrap(~RCP,ncol=1)
dev.off()


range(tapply(NovSull$Frec,
list(NovSull$Biome,
     NovSull$Mech,
     NovSull$RCP),
sum)[,"NovClim/FastDis/OrthDiv",])

apply(tapply(NovSull$Frec,
             list(NovSull$Biome,
                  NovSull$Mech,
                  NovSull$RCP),
             sum)[-1,"NovClim/FastDis/OrthDiv",],
      2,
      median)


apply(tapply(NovSull$Frec,
             list(NovSull$Biome,
                  NovSull$Mech,
                  NovSull$RCP),
             sum)[-1,"NovClim/FastDis",],
      2,
      median)

