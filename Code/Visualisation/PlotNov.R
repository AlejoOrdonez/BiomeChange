rm(list=ls());gc()
require(terra)
library(tidyverse)
library(tidyterra)
require(maptools)
data(wrld_simpl)
wrld_simpl2 <- vect(wrld_simpl)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange/Results/Displacement_Divergence/Seasonal")