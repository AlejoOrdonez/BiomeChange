rm(list=ls());gc()
require(terra)
library(tidyterra)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sp")
wrld_simpl2 <- vect(world)
wrld_simpl2 <- project(wrld_simpl2,"+proj=eck4")
rm(world);gc()

# Estimate leverage of each model for Climate novelty Average
Climate_novelty_Model_Leverage <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                                         function(RCP){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[1])
                                           Climate_novelty_Model_Leverage <- lapply(seq(2100,2300,by=50),
                                                                                    function(YearUse){#(YearUse <- seq(2100,2300,by=50)[2])
                                                                                      # Define the models run for a given scenario
                                                                                      Models_run <- unique(sapply(strsplit(dir(paste0("./Results/Novelty/By_Model/",RCP)),"_"),
                                                                                                                  function(x){x[2]}))
                                                                                      # Load the raster of Climate novelty for each model
                                                                                      Climate_novelty_bymodel_List<- lapply(Models_run,
                                                                                                                            function(ModelUse){#(ModelUse <- Models_run[1])
                                                                                                                              # Analog threhold Load
                                                                                                                              Tresh <- readRDS(paste0("./Results/Novelty/By_Model/Tresholds/TreshSumm_",ModelUse,".rds"))$MDSummTresh
                                                                                                                              # Load the Disimilarity
                                                                                                                              MDRast <- rast(paste0("./Results/Novelty/By_Model/",RCP,"/MDminSumm_",ModelUse,"_",RCP,"_",YearUse,".tif"),
                                                                                                                                   lyrs = 1)
                                                                                                                              NoAnalogue <- (MDRast>Tresh)*1
                                                                                                                              return(NoAnalogue)
                                                                                                                            })
                                                                                      # Turn the list into a raster
                                                                                      Climate_novelty_bymodel_rast <- do.call("c",Climate_novelty_bymodel_List)
                                                                                      names(Climate_novelty_bymodel_rast) <- Models_run
                                                                                      AggModel <- 100*(sum(Climate_novelty_bymodel_rast)/length(Models_run))
                                                                                      return(AggModel)})
                                           Climate_novelty_Model_Leverage <- do.call("c",Climate_novelty_Model_Leverage)
                                           names(Climate_novelty_Model_Leverage) <- seq(2100,2300,by=50)
                                           return(Climate_novelty_Model_Leverage)
                                         })
Climate_novelty_Model_Leverage <- do.call("c",Climate_novelty_Model_Leverage)
names(Climate_novelty_Model_Leverage) <- paste0(rep(c("RCP26", "RCP45", "RCP60", "RCP85"),each=5),"-",names(Climate_novelty_Model_Leverage) )
# Plot the leverage of each model on the cross model summary
pdf("./Results/PDF/Sup_Figure_10.pdf",
    width=15, height =10)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = Climate_novelty_Model_Leverage) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #color scheme 
                        na.value = NA,#Do not map NA
                        breaks = seq(0,100,by=10), # Legend breaks
                        name = "Agrement in Novelty",
                        limit = c(0,100))+ # Legend Title
  facet_wrap(~lyr, ncol=5) +
  theme(legend.position="bottom",legend.key.width = unit(4, 'cm')) +
  geom_spatvector(fill = NA,linewidth = 0.2) # Add the vector of the world
#scale_x_discrete(guide = guide_axis(n.dodge = 2))
dev.off()

# Estimate leverage of each model for Displacement novelty average
Climate_Disp_Model_Leverage <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                                        function(RCP){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[1])
                                          Climate_Disp_Model_Leverage <- lapply(seq(2100,2300,by=50),
                                                                                  function(YearUse){#(YearUse <- seq(2100,2300,by=50)[2])
                                                                                    # Define the models run for a given scenario
                                                                                    Models_run <- unique(sapply(strsplit(dir(paste0("./Results/Displacement_Divergence/By_Model/",RCP)),"_"),
                                                                                                                function(x){x[3]}))
                                                                                    # Load the raster of Climate novelty for each model
                                                                                    Climate_Disp_bymodel_List<- lapply(Models_run,
                                                                                                                         function(ModelUse){#(ModelUse<-Models_run[1])
                                                                                                                           MDRast <- rast(paste0("./Results/Displacement_Divergence/By_Model/",RCP,"/Displace_Diverg_",ModelUse,"_",RCP,"_",YearUse,".tif"),
                                                                                                                                lyrs = 1)
                                                                                                                           NoAnalogue <- (MDRast>1)*1
                                                                                                                           return(NoAnalogue)
                                                                                                                         })
                                                                                    # Turn the list into a raster
                                                                                    Climate_Disp_bymodel_rast <- do.call("c",Climate_Disp_bymodel_List)
                                                                                    names(Climate_Disp_bymodel_rast) <- Models_run
                                                                                    AggModel <- 100*(sum(Climate_Disp_bymodel_rast)/length(Models_run))
                                                                                    return(AggModel)})
                                          Climate_Disp_Model_Leverage <- do.call("c",Climate_Disp_Model_Leverage)
                                          names(Climate_Disp_Model_Leverage) <- seq(2100,2300,by=50)
                                          return(Climate_Disp_Model_Leverage)
                                        })
Climate_Disp_Model_Leverage <- do.call("c",Climate_Disp_Model_Leverage)
names(Climate_Disp_Model_Leverage) <- paste0(rep(c("RCP26", "RCP45", "RCP60", "RCP85"),each=5),"-",names(Climate_Disp_Model_Leverage) )
# Plot the leverage of each model on the cross model summary
pdf("./Results/PDF/Sup_Figure_11.pdf",
    width=15, height =10)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = Climate_Disp_Model_Leverage) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #color scheme 
                        na.value = NA,#Do not map NA
                        breaks = seq(0,100,by=10), # Legend breaks
                        name = "Agrement in Novelty",
                        limit = c(0,100))+ # Legend Title
  facet_wrap(~lyr, ncol=5) +
  theme(legend.position="bottom",legend.key.width = unit(4, 'cm')) +
  geom_spatvector(fill = NA,linewidth = 0.2) # Add the vector of the world
#scale_x_discrete(guide = guide_axis(n.dodge = 2))
dev.off()


# Estimate leverage of each model for Divergence novelty Average
Climate_Diverg_Model_Leverage <- lapply(c("RCP26", "RCP45", "RCP60", "RCP85"),
                                        function(RCP){#(RCP <- c("RCP26", "RCP45", "RCP60", "RCP85")[1])
                                          Climate_Diverg_Model_Leverage <- lapply(seq(2100,2300,by=50),
                                                                                  function(YearUse){#(YearUse <- seq(2100,2300,by=50)[2])
                                                                                    # Define the models run for a given scenario
                                                                                    Models_run <- unique(sapply(strsplit(dir(paste0("./Results/Displacement_Divergence/By_Model/",RCP)),"_"),
                                                                                                                function(x){x[3]}))
                                                                                    # Load the raster of Climate novelty for each model
                                                                                    Climate_Diverg_bymodel_List<- lapply(Models_run,
                                                                                                                         function(ModelUse){#(ModelUse<-Models_run[1])
                                                                                                                           MDRast <-rast(paste0("./Results/Displacement_Divergence/By_Model/",RCP,"/Displace_Diverg_",ModelUse,"_",RCP,"_",YearUse,".tif"),
                                                                                                                                lyrs = 2)
                                                                                                                           NoAnalogue <- (MDRast > 60 & MDRast < 120)*1
                                                                                                                           return(NoAnalogue)
                                                                                                                         })
                                                                                    # Turn the list into a raster
                                                                                    Climate_Diverg_bymodel_rast <- do.call("c",Climate_Diverg_bymodel_List)
                                                                                    names(Climate_Diverg_bymodel_rast) <- Models_run
                                                                                    AggModel <- 100*(sum(Climate_Diverg_bymodel_rast)/length(Models_run))
                                                                                    return(AggModel)})
                                          Climate_Diverg_Model_Leverage <- do.call("c",Climate_Diverg_Model_Leverage)
                                          names(Climate_Diverg_Model_Leverage) <- seq(2100,2300,by=50)
                                          return(Climate_Diverg_Model_Leverage)
                                        })
Climate_Diverg_Model_Leverage <- do.call("c",Climate_Diverg_Model_Leverage)
names(Climate_Diverg_Model_Leverage) <- paste0(rep(c("RCP26", "RCP45", "RCP60", "RCP85"),each=5),"-",names(Climate_Diverg_Model_Leverage) )
# Plot the leverage of each model on the cross model summary
pdf("./Results/PDF/Sup_Figure_12.pdf",
    width=15, height =10)
ggplot(wrld_simpl2) + # add the vector of the world
  geom_spatraster(data = Climate_Diverg_Model_Leverage) + # Map the Displacement
  # Setup. plot of a continuous raster
  scale_fill_whitebox_c(palette = "muted", #color scheme 
                        na.value = NA,#Do not map NA
                        breaks = seq(0,100,by=10), # Legend breaks
                        name = "Agrement in Novelty",
                        limit = c(0,100))+ # Legend Title
  facet_wrap(~lyr, ncol=5) +
  theme(legend.position="bottom",legend.key.width = unit(4, 'cm')) +
  geom_spatvector(fill = NA,linewidth = 0.2) # Add the vector of the world
#scale_x_discrete(guide = guide_axis(n.dodge = 2))
dev.off()