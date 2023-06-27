rm(list=ls());gc()
library(terra)
setwd("~/Downloads/AllModels")


# Analysis for the RCP-2.6
# Load the analogue treshold
Tresh <- readRDS("./RCP26/BIOCLIM/AllModels_RCP26_TreshSumm.rds")

# Load the Novelty metrics for 2100
RCP26.2100 <- rast("./RCP26/BIOCLIM/AllModels_RCP26_2099_MDminSumm.tif")
# Plot novel areas by 2100
plot(RCP26.2100[[1]]>Tresh$MDSummTresh)
# How many cells are novel
sum(values(RCP26.2100[[1]]>Tresh$MDSummTresh,na.rm=T)[,1])
# % of the arteh to be novel
sum(values(RCP26.2100[[1]]>Tresh$MDSummTresh,na.rm=T)[,1])/dim(values(RCP26.2100[[1]],na.rm=T))[1]

# Load the Novelty metrics for 2300
RCP26.2300 <- rast("./RCP26/BIOCLIM/AllModels_RCP26_2299_MDminSumm.tif")
# Plot novel areas by 2100
plot(RCP26.2300[[1]]>Tresh$MDSummTresh)
# How many cells are novel
sum(values(RCP26.2300[[1]]>Tresh$MDSummTresh,na.rm=T)[,1])
# % of the arteh to be novel
sum(values(RCP26.2300[[1]]>Tresh$MDSummTresh,na.rm=T)[,1])/dim(values(RCP26.2100[[1]],na.rm=T))[1]


# Analysis for the RCP-8.5
# Load the analogue treshold
Tresh <- readRDS("./RCP85/BIOCLIM/AllModels_RCP85_TreshSumm.rds")

# Load the Novelty metrics for 2100
RCP85.2100 <- rast("./RCP85/BIOCLIM/AllModels_RCP85_2099_MDminSumm.tif")
# Plot novel areas by 2100
plot(RCP85.2100[[1]]>Tresh$MDSummTresh)
# How many cells are novel
sum(values(RCP85.2100[[1]]>Tresh$MDSummTresh,na.rm=T)[,1])
# % of the arteh to be novel
sum(values(RCP85.2100[[1]]>Tresh$MDSummTresh,na.rm=T)[,1])/dim(values(RCP85.2100[[1]],na.rm=T))[1]

# Analogue based velocity
plot((RCP85.2100[[3]]/100))
Dist.2100 <- as.numeric(values(RCP85.2100[[3]],na.rm=T))
Dist.2100 <- Dist.2100[Dist.2100!=0]
# Avg Speed - anlaogue based
10^mean(log10(Dist.2100/100))
hist(log10(Dist.2100/100))
abline(v=mean(log10(Dist.2100/100)), col="red")




# Load the Novelty metrics for 2300
RCP85.2300 <- rast("./RCP85/BIOCLIM/AllModels_RCP85_2299_MDminSumm.tif")
# Plot novel areas by 2100
plot(RCP85.2300[[1]]>Tresh$MDSummTresh)
# How many cells are novel
sum(values(RCP85.2300[[1]]>Tresh$MDSummTresh,na.rm=T)[,1])
# % of the arteh to be novel
sum(values(RCP85.2300[[1]]>Tresh$MDSummTresh,na.rm=T)[,1])/dim(values(RCP85.2100[[1]],na.rm=T))[1]

# Analogue based velocity
plot(log10(RCP85.2300[[3]]/100))
Dist.2300 <- as.numeric(values(RCP85.2300[[3]],na.rm=T))
Dist.2300 <- Dist.2300[Dist.2300!=0]
# Avg Speed - anlaogue based
10^mean(log10(Dist.2300/300))
hist(log10(Dist.2300/300))
abline(v=mean(log10(Dist.2300/300)), col="red")
