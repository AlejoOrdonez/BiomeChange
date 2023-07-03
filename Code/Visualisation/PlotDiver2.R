rm(list=ls());gc()
require(terra)
library(tidyverse)
library(tidyterra)

# Load a simpified map of the world
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


# Plot the divergence for each of the evaluated time periods (2049,2099,2149,2199,2249,2299) for all RCPs
# Displacement estimated using an ARM model
Model <- c("BIOCLIM","BIOCLIM2","koeppen_geiger","Seasonal")[4]#(Model <- c("BIOCLIM","BIOCLIM2","koeppen_geiger","Seasonal")[1])

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

# Map
pdf("/Volumes/MacPro 2013 Backup/BiomeChange/Results2/PDF/Fig6.pdf",
    width = 3, height=9)#width = 10, height=5)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = DivergenceSumm) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = c(0,60,120,180), # Legend breaks
                        #labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Divergence angle", # Legend Title
                        limit = c(0,180)) +
  facet_wrap(~lyr,ncol = 1,dir= "v") +
  theme(legend.position="bottom") +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = Model) # Fig title
dev.off()


# mean divergence of median divergence accords periods
lapply(Divergence,
       function(x){#(x <- Divergence[[1]])
         ValuesTmp <- values(x,na.rm=T)
         mean(apply(ValuesTmp,2,median))})
# mean divergence of mean divergence accords periods
lapply(Divergence,
       function(x){#(x <- Divergence[[1]])
         ValuesTmp <- values(x,na.rm=T)
         mean(ValuesTmp)})


# Proportion of cells with Congruent (0-60), orthogonal (60-120), and opposing (120-180) divergences
DivergenceTbl1 <- lapply(Divergence,
                         function(x){#(x <- Divergence[[1]])
                           ValuesTmp <- values(x,na.rm=T)
                           DivTmp <- apply(ValuesTmp,
                                           2,
                                           function(y){#{y <- ValuesTmp[,y]}
                                             y <-factor(ifelse(y<60,1,ifelse(y<120,2,3)),
                                                        levels = 1:3)
                                             y <- as.numeric(round(100*(table(y)/length(y)),2))
                                             names(y) <- c("congruent","orthogonal","opposing")
                                             return(y)
                                           })
                           return(DivTmp)
                         })

DivergenceTbl2 <- data.frame(RCP = factor(rep(names(DivergenceTbl1),each=length(DivergenceTbl1[[1]]))),
                             Year = factor(rep(colnames(DivergenceTbl1[[1]]),each=3)),
                             Synd = factor(rep(rownames(DivergenceTbl1[[1]]),6)),
                             Freq = do.call("c",(lapply(DivergenceTbl1,function(x){c(x)}))))


ggplot(DivergenceTbl2, aes(fill=Synd, y=Freq, x=Year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "RdYlBu") +
  ylab("Proportion Land under a divergecne syndrome") +
  xlab("Time period") +
  facet_wrap(~RCP) +
  ggtitle(Model)



# Summary of divergences per Biome
DivSumByBiome <- lapply(DivergenceSumm,
                        function(x){#(x <- DivergenceSumm[[1]])
                          a <- terra::extract(x, BiomeShp, fun = mean,na.rm=T)
                          a <- tapply(a[,2],BiomeShp$ECO_ID,mean,na.rm=T)
                          a.1 <- tapply(a,
                                        BiomeShp$BIOME[match(names(a),BiomeShp$ECO_ID)],
                                        mean,na.rm=T)
                          names(a.1) <- BiomeNames$Name
                          return(a.1)
                        })

range(sapply(DivSumByBiome,range))

DivSumByBiomeTbl <- data.frame(RCP = rep(names(DivergenceSumm), each = length(DivSumByBiome[[1]])),
                               Biome = (names(unlist(DivSumByBiome))),
                               Div = as.numeric(unlist(DivSumByBiome)))
ggplot(DivSumByBiomeTbl, aes(y=Div, x=Biome)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("Divergence") +
  xlab("Biome") +
  scale_x_discrete(limits= BiomeNames$Name) +
  ylim(0,180) +
  ggtitle(Model) +
  geom_hline(yintercept=c(60,90,120), linetype="dashed") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  facet_wrap(~RCP)



# Summary of divergences per Biome
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

pdf("/Volumes/MacPro 2013 Backup/BiomeChange/Results2/PDF/Fig7.pdf",
    width = 4, height=10)
ggplot(DivSumByBiome2, aes(fill=Sind, y=Freq, x=Biome)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "RdYlBu") +
  ylab("Proportion Land under a divergecne syndrome") +
  xlab("Time period") +
  facet_wrap(~RCP,ncol=1) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  ggtitle(Model)
dev.off()