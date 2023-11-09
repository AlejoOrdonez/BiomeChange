rm(list=ls());gc()
require(terra)
require(snowfall)
require(analogue)
require(tidyverse)

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
                                                                                                                            function(ModelUse){
                                                                                                                              rast(paste0("./Results/Novelty/By_Model/",RCP,"/MDminSumm_",ModelUse,"_",RCP,"_",YearUse,".tif"),
                                                                                                                                   lyrs = 1)
                                                                                                                            })
                                                                                      # Turn the list into a raster
                                                                                      Climate_novelty_bymodel_rast <- do.call("c",Climate_novelty_bymodel_List)
                                                                                      names(Climate_novelty_bymodel_rast) <- Models_run
                                                                                      
                                                                                      # estimate the leverage of each model in the Novelty esttimates
                                                                                      Climate_novelty_Model_Leverage <- lapply(1:dim(Climate_novelty_bymodel_rast)[3],
                                                                                                                               function(x){
                                                                                                                                 
                                                                                                                                 #Define the Median of all the Models
                                                                                                                                 Median_models <- app(Climate_novelty_bymodel_rast,median)
                                                                                                                                 
                                                                                                                                 #Define the Median with out the ith- Model
                                                                                                                                 Reduced_model <- app(Climate_novelty_bymodel_rast[[-x]],median)
                                                                                                                                 
                                                                                                                                 # Estimate the devance betwen the fuall and one model reeuded estimate of novelty
                                                                                                                                 data.frame(Models = Models_run[x],
                                                                                                                                            median = log10(values(Reduced_model-Median_models,df=T,na.rm=T)+1))
                                                                                                                               })
                                                                                      
                                                                                      Climate_novelty_Model_Leverage <- do.call("rbind",Climate_novelty_Model_Leverage)
                                                                                      Climate_novelty_Model_Leverage$Year <- YearUse
                                                                                      Climate_novelty_Model_Leverage$RCP <- RCP
                                                                                      return(Climate_novelty_Model_Leverage)
                                                                                    })
                                           Climate_novelty_Model_Leverage <- do.call("rbind",Climate_novelty_Model_Leverage)
                                           return(Climate_novelty_Model_Leverage)
                                         })
Climate_novelty_Model_Leverage <- do.call("rbind",Climate_novelty_Model_Leverage)
Climate_novelty_Model_Leverage$Year <- factor(Climate_novelty_Model_Leverage$Year)

# Plot the leverage of each model on the cross model summary
pdf("./Results/PDF/Sup_Figure_7.pdf",
    width=15, height =10)
ggplot(Climate_novelty_Model_Leverage) + aes(x=median,fill=Models) + geom_density(alpha=0.3) + facet_wrap(~RCP*Year) +xlim(-0.025,0.025) + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
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
                                                                                                                           rast(paste0("./Results/Displacement_Divergence/By_Model/",RCP,"/Displace_Diverg_",ModelUse,"_",RCP,"_",YearUse,".tif"),
                                                                                                                                lyrs = 1)
                                                                                                                         })
                                                                                    # Turn the list into a raster
                                                                                    Climate_Disp_bymodel_rast <- do.call("c",Climate_Disp_bymodel_List)
                                                                                    names(Climate_Disp_bymodel_rast) <- Models_run
                                                                                    
                                                                                    # estimate the leverage of each model in the Novelty esttimates
                                                                                    Climate_Disp_Model_Leverage <- lapply(1:dim(Climate_Disp_bymodel_rast)[3],
                                                                                                                            function(x){
                                                                                                                              
                                                                                                                              #Define the Median of all the Models
                                                                                                                              Median_models <- app(Climate_Disp_bymodel_rast,median)
                                                                                                                              
                                                                                                                              #Define the Median with out the ith- Model
                                                                                                                              Reduced_model <- app(Climate_Disp_bymodel_rast[[-x]],median)
                                                                                                                              
                                                                                                                              # Estimate the devance betwen the fuall and one model reeuded estimate of novelty
                                                                                                                              data.frame(Models = Models_run[x],
                                                                                                                                         median = values(Reduced_model-Median_models,df=T,na.rm=T))
                                                                                                                            })
                                                                                    
                                                                                    Climate_Disp_Model_Leverage <- do.call("rbind",Climate_Disp_Model_Leverage)
                                                                                    Climate_Disp_Model_Leverage$Year <- YearUse
                                                                                    Climate_Disp_Model_Leverage$RCP <- RCP
                                                                                    return(Climate_Disp_Model_Leverage)
                                                                                  })
                                          Climate_Disp_Model_Leverage <- do.call("rbind",Climate_Disp_Model_Leverage)
                                          return(Climate_Disp_Model_Leverage)
                                        })
Climate_Disp_Model_Leverage <- do.call("rbind",Climate_Disp_Model_Leverage)
Climate_Disp_Model_Leverage$Year <- factor(Climate_Disp_Model_Leverage$Year)

# Plot the leverage of each model on the cross model summary
pdf("./Results/PDF/Sup_Figure_8.pdf",
    width=15, height =10)
ggplot(Climate_Disp_Model_Leverage) + aes(x=median,fill=Models) + geom_density(alpha=0.3) + facet_wrap(~RCP*Year) +xlim(-0.2,0.2) + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
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
                                                                                                                           rast(paste0("./Results/Displacement_Divergence/By_Model/",RCP,"/Displace_Diverg_",ModelUse,"_",RCP,"_",YearUse,".tif"),
                                                                                                                                lyrs = 2)
                                                                                                                         })
                                                                                    # Turn the list into a raster
                                                                                    Climate_Diverg_bymodel_rast <- do.call("c",Climate_Diverg_bymodel_List)
                                                                                    names(Climate_Diverg_bymodel_rast) <- Models_run
                                                                                    
                                                                                    # estimate the leverage of each model in the Novelty esttimates
                                                                                    Climate_Diverg_Model_Leverage <- lapply(1:dim(Climate_Diverg_bymodel_rast)[3],
                                                                                                                            function(x){
                                                                                                                              
                                                                                                                              #Define the Median of all the Models
                                                                                                                              Median_models <- app(Climate_Diverg_bymodel_rast,median)
                                                                                                                              
                                                                                                                              #Define the Median with out the ith- Model
                                                                                                                              Reduced_model <- app(Climate_Diverg_bymodel_rast[[-x]],median)
                                                                                                                              
                                                                                                                              # Estimate the devance betwen the fuall and one model reeuded estimate of novelty
                                                                                                                              data.frame(Models = Models_run[x],
                                                                                                                                         median = values(Reduced_model-Median_models,df=T,na.rm=T))
                                                                                                                            })
                                                                                    
                                                                                    Climate_Diverg_Model_Leverage <- do.call("rbind",Climate_Diverg_Model_Leverage)
                                                                                    Climate_Diverg_Model_Leverage$Year <- YearUse
                                                                                    Climate_Diverg_Model_Leverage$RCP <- RCP
                                                                                    return(Climate_Diverg_Model_Leverage)
                                                                                  })
                                          Climate_Diverg_Model_Leverage <- do.call("rbind",Climate_Diverg_Model_Leverage)
                                          return(Climate_Diverg_Model_Leverage)
                                        })
Climate_Diverg_Model_Leverage <- do.call("rbind",Climate_Diverg_Model_Leverage)
Climate_Diverg_Model_Leverage$Year <- factor(Climate_Diverg_Model_Leverage$Year)

# Plot the leverage of each model on the cross model summary
pdf("./Results/PDF/Sup_Figure_9.pdf",
    width=15, height =10)
ggplot(Climate_Diverg_Model_Leverage) + aes(x=median,fill=Models) + geom_density(alpha=0.3) + facet_wrap(~RCP*Year) +xlim(-3.5,3.5) + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
#scale_x_discrete(guide = guide_axis(n.dodge = 2))
dev.off()
