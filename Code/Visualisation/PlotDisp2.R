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


# Plot the Displacement for each of the evaluated time periods (2049,2099,2149,2199,2249,2299) for all RCPs
# Displacement estimated using an ARM model
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


# Plot the Displacement for RCP 8.5
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = log10(DisplARM[["RCP85"]])) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = -2:2, # Legend breaks
                        labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Speed\n[km/yrs]", # Legend Title
                        limit = c(-2,2)) +
  facet_wrap(~lyr) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = "Displacement RCP 8.5") # Fig title


# Plot the Displacement for RCP 8.5 vs RCP 2.6

DisplARM_26vs85 <- c(DisplARM[["RCP85"]],DisplARM[["RCP26"]])
names(DisplARM_26vs85) <- paste0(rep(c("RCP 8.5","RCP 2.6"),each=6),"-",
                                 names(DisplARM_26vs85))

# Mean difference in displacement per cell between RCP 2.6 and RCP8.5
mean(c(mean(values(DisplARM_26vs85[[1]]-DisplARM_26vs85[[7]]),na.rm=T),
       mean(values(DisplARM_26vs85[[2]]-DisplARM_26vs85[[8]]),na.rm=T),
       mean(values(DisplARM_26vs85[[3]]-DisplARM_26vs85[[9]]),na.rm=T),
       mean(values(DisplARM_26vs85[[4]]-DisplARM_26vs85[[10]]),na.rm=T),
       mean(values(DisplARM_26vs85[[5]]-DisplARM_26vs85[[11]]),na.rm=T),
       mean(values(DisplARM_26vs85[[6]]-DisplARM_26vs85[[12]]),na.rm=T)))


# Plot the projected RCP 2.6 and RCP8.5 displacment 
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = log10(DisplARM_26vs85)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = -2:2, # Legend breaks
                        labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Speed\n[km/yrs]", # Legend Title
                        limit = c(-2,2)) +
  facet_wrap(~lyr,ncol = 2,dir= "v") +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = "Displacement RCP 8.5 vs RCP 2.6") # Fig title


# Plot the Displacement for the 2049 to 2099 period
Displ_2049to2099 <- c(DisplARM[["RCP26"]][["2049 to 2099"]],
                      DisplARM[["RCP45"]][["2049 to 2099"]],
                      DisplARM[["RCP60"]][["2049 to 2099"]],
                      DisplARM[["RCP85"]][["2049 to 2099"]])
names(Displ_2049to2099) <- c("RCP26", "RCP45", "RCP60", "RCP85")
# Plot the Velocities 2049
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = log10(Displ_2049to2099)) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #colour scheme 
                        na.value = NA,#Do not map NA
                        breaks = -2:2, # Legend breaks
                        labels = c("<0.01","0.1","1","10",">100"), # Legend Labels
                        name = "Speed\n[km/yrs]", # Legend Title
                        limit = c(-2,2)) +
  facet_wrap(~lyr) +
  geom_spatvector(fill = NA) + # Add the vector of the world
  labs(title = "Displacement 2049 to 2099") # Fig title


## Summary of Displacements globaly
### Calculate the geometric mean of the displacement per Yr/RCP
DisARMSumm <- lapply(DisplARM,
                     function(x){
                       out<- data.frame(Year = factor(names(x)),
                                  DisplArm = 10^apply(log10(values(x,na.rm=T)),
                                                      2,
                                                      mean))
                     })
### Compile geometric means of the displacement per Yr/RCP as a dataframe
DisARMSumm2 <- do.call("rbind",DisARMSumm)
DisARMSumm2<- data.frame(RCP = rep(c("RCP26", "RCP45", "RCP60", "RCP85"),
                                    each=dim(DisARMSumm[[1]])[1]),
                         DisARMSumm2)
# Plot the change in 50yrs velocities per RCP
ggplot(DisARMSumm2,
      aes(x=Year, y=DisplArm, group=RCP,color = RCP)) + 
  ylab("Displaecment\n[km/year]") +
  xlab("Time Period")+
  geom_line() +
  geom_point(shape = 19) +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values=c("darkgreen","purple","gold","red")) +
  geom_hline(yintercept=1, linetype="dashed", color = "red")



# Decrease in displacment per year
summary(lm(DisARMSumm2[DisARMSumm2$RCP=="RCP26","DisplArm"][1:3]~seq(2049,2299,by=50)[1:3]))
summary(lm(DisARMSumm2[DisARMSumm2$RCP=="RCP45","DisplArm"][1:3]~seq(2049,2299,by=50)[1:3]))

summary(lm(DisARMSumm2[DisARMSumm2$RCP=="RCP60","DisplArm"][2:4]~seq(2049,2299,by=50)[2:4]))
summary(lm(DisARMSumm2[DisARMSumm2$RCP=="RCP85","DisplArm"][2:6]~seq(2049,2299,by=50)[2:6]))




### Calculate the geometric mean of the displacement per Biome and Yr/RCP
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


DisARMBiomeSumm <- data.frame (RCP = rep(names(DisplARM), each = dim(DisARMBiomeSumm[[1]])[1]),
                               do.call("rbind",DisARMBiomeSumm))
#DisARMBiomeSumm$DisplArm <- log10(DisARMBiomeSumm$DisplArm)

# plot the Rate of Displacement 
ggplot(DisARMBiomeSumm,
       aes(x=Year, y=DisplArm, group=Biome ,color = Biome)) + 
  ylab("Displaecment\n[km/year]") +
  xlab("Time Period") +
  geom_line() +
  geom_point(shape = 19) +
  scale_y_continuous(trans='log10') + 
  #scale_color_manual(values=c("darkgreen","purple","gold","red")) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  facet_wrap(facets = "RCP")



