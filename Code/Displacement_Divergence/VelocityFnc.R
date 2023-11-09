### FUNCTIONS VELOCITY - DISPLACEMENT - DIVERGENCE
### BY: Alejandro Ordonez
### Date: 10th Feb 2023
#-------------------------------------------------------------------------------
## TempGradFnc: Function to estimate temporal gradients. For this, you need to specify:
####  The vector input a row in a raster RastIn,
####  The temporal spacing in years between raster layers [TimeStep]
####  The Method Used to estimate the Change over time [method "ARM","glm","lm","AbsDif","Anomaly"]
####  The regression family to use if method is glm [FamilyUse]
TempGradFnc <- function(x,
                        method = "ARM",
                        FamilyUse = "binomial",
                        TimeStep = 1){
  if(!method%in%c("ARM","glm","lm","AbsDif","Anomaly1","Anomaly2")){
    stop("method has to be: ARM, glm, lm, AbsDif, Anomaly1, or Anomaly2")
  }
  if (sum(x,na.rm=T)!=0 & sum(x>0,na.rm=T)>3 & length(unique(x)) > 2){# & length(unique(x))>3){
    tmpDtaFrm <- data.frame(prop = x[1:length(x)],
                            Time = c(1:length(x)))
    
    if(method == "ARM"){ # Slope of a regression with temporal autocorrelation structure --> AR1 model
      require(nlme)
      TimMod <- gls(prop~Time, data = tmpDtaFrm,
                    correlation = corARMA(p=1), 
                    method ="ML")
      if(!"TimMod"%in%ls()){
      TimMod <- gls(prop~Time, data = tmpDtaFrm, 
                    method ="ML")}
      Out <- coef(summary(TimMod))["Time",c("Value")]/(TimeStep)#c("Value","Std.Error","p-value")
    }
    if(method == "glm"){ # Slope of a GLM regression model of a given Family
      TimMod <- glm(prop~Time,
                    data = tmpDtaFrm,
                    family = FamilyUse) 
      Out <- coef(summary(TimMod))["Time",c("Estimate")]/(TimeStep)#c("Estimate","Std. Error","Pr(>|z|)")
    }
    if(method == "lm"){ # Slope of a GLS regression model
      TimMod <- lm(prop~Time, data = tmpDtaFrm)
      Out <- coef(summary(TimMod))["Time",c("Estimate")]/(TimeStep)#c("Estimate","Std. Error","Pr(>|t|)")
    }
    if(method == "AbsDif"){ # Total difference between consecutive time periods divided by the number of years
      TimMod <- tmpDtaFrm$prop[-1] - tmpDtaFrm$prop[-dim(tmpDtaFrm)[1]]
      Out <- sum(abs(TimMod))/((TimeStep*(length(x)-1)))
    }
    if(method == "Anomaly1"){ # median difference between consecutive time periods
      TimMod <- tmpDtaFrm$prop[-1] - tmpDtaFrm$prop[-dim(tmpDtaFrm)[1]]
      Out <- median(TimMod)/(TimeStep)
    }
    if(method == "Anomaly2"){ # Difference between start and end period divided by the number of years
      Out <- (tmpDtaFrm$prop[1] - tmpDtaFrm$prop[dim(tmpDtaFrm)[1]])/((TimeStep*(length(x)-1)))
    }
  }
  else{
    Out <- 0#c(0,0,1)
  }
  return(Out)
}
#-------------------------------------------------------------------------------
## SpatHetFnc: Function to estimate the spatial gradients magnitude using a
##             using the maximum average technique [Burrough & McDonnell 1998].
####  The Raster input [RastIn],
SpatHetFnc <- function(RastIn,Dist){
  if(dim(RastIn)[3]!=1){
    RastIn <- mean(RastIn)
    warning("Input raster has more than one layer - and averaged raster is used")
  }
  # West-east gradients
  EstWestChngEvrGrn <- focal(RastIn,# Input Raster
                             w = 3, # Neighborhood matrix
                             fun = function(x){mean(c(x[2]-x[1],x[3]-x[2],
                                                      x[5]-x[4],x[6]-x[5],
                                                      x[8]-x[7],x[9]-x[8]),
                                                    na.rm=TRUE)/Dist
                             })
  
  # Poleward gradients - Norther hemisphere (negative change means equatorial movement)
  RastInNorth <- crop(RastIn,
                      ext(as.numeric(c(ext(RastIn)[c(1,2)],
                                       0,
                                       ext(RastIn)[4]))))
  
  NrthSthChngEvrGrn1 <- focal(RastIn,# Input Raster
                              w = 3, # Neighborhood matrix
                              fun = function(x){mean(c(x[1]-x[4],x[4]-x[7],
                                                       x[2]-x[5],x[5]-x[8],
                                                       x[3]-x[6],x[6]-x[9]),
                                                     na.rm=TRUE)/Dist
                              })
  # Poleward gradients - South hemisphere  (negative change means equatorial movement)
  RastInSouth <- crop(RastIn,
                      ext(as.numeric(c(ext(RastIn)[c(1,2,3)],0))))
  
  NrthSthChngEvrGrn2 <- focal(RastIn,# Input Raster
                              w = 3, # Neighborhood matrix
                              fun = function(x){mean(c(x[4]-x[1],x[7]-x[4],
                                                       x[5]-x[2],x[8]-x[5],
                                                       x[6]-x[3],x[9]-x[6]),
                                                     na.rm=TRUE)/Dist
                              })
  # Mosaic the North-South gradients 
  NrthSthChngEvrGrn <- mosaic(NrthSthChngEvrGrn1,NrthSthChngEvrGrn2)
  
  # Vector sum of the N-S and W-E gradients --> the magnitude of the spatial gradient
  SpacHetEvGrnRast <- sqrt((NrthSthChngEvrGrn^2)+(EstWestChngEvrGrn^2))
  return(SpacHetEvGrnRast)
}
#-------------------------------------------------------------------------------
## BearingFnc: Function to estimate the spatial gradients bearing using a
##             using the maximum average technique [Burrough & McDonnell 1998].
##             The angle is define in a poleward direction, so 90-D point to the
##             pole. Also 361 means that there is no spatial gradient (no change
##            in the x or y direction.)
####  The Raster input [RastIn],
BearingFnc <- function(RastIn){
  # Check for right dimensionality
  if(dim(RastIn)[3]==1){
    stop("Input raster needs at least to layers [strat and end period]")
  } #RastIn <- MapsPer.GF
  # Time integrated bearing The direction of change would change if the temporal tendency is to a reduction in the variable of interest
  # Load the temporal trend  
  ### Positive value means that the past is more suitable than the present so the direction is out of the cell
  ### Negative value means that the present is more suitable than the past so the direction is into of the cell
  TempHetRast <- RastIn[[1]] - RastIn[[dim(RastIn)[3]]]
  # Define the Baseline raster
  RastUse <- RastIn[[1]]#RastUse <- mean(MapsPer.GF)#
  # West-east gradients
  EstWestChng <- focal(RastUse,# Input Raster
                       w = 3, # Neighborhood matrix
                       fun = function(x){mean(c(x[2]-x[1],x[3]-x[2],
                                                x[5]-x[4],x[6]-x[5],
                                                x[8]-x[7],x[9]-x[8]),
                                              na.rm=TRUE)/47
                       })
  # North-South gradients - Norther hemisphere (negative change means equatorial movement)
  RastIn.North <- crop(RastUse,
                       ext(as.numeric(c(ext(RastIn)[c(1,2)],
                                        0,
                                        ext(RastIn)[4]))))
  NrthSthChng1 <- focal(RastIn.North,# Input Raster
                        w = 3, # Neighborhood matrix
                        fun = function(x){mean(c(x[1]-x[4],x[4]-x[7],
                                                 x[2]-x[5],x[5]-x[8],
                                                 x[3]-x[6],x[6]-x[9]),
                                               na.rm=TRUE)/65.9
                        })
  # North-South gradients - South hemisphere  (negative change means equatorial movement)
  RastIn.South <- crop(RastUse,
                       ext(as.numeric(c(ext(RastIn)[c(1,2,3)],0))))
  NrthSthChng2 <- focal(RastIn.South,# Input Raster
                        w = 3, # Neighborhood matrix
                        fun = function(x){mean(c(x[4]-x[1],x[7]-x[4],
                                                 x[5]-x[2],x[8]-x[5],
                                                 x[6]-x[3],x[9]-x[6]),
                                               na.rm=TRUE)/65.9
                        })
  # Mosaic the North-South gradients 
  NrthSthChng<- mosaic(NrthSthChng1,NrthSthChng2)
# Estimating the bearing of the velocity vector based on initial conditions
  BearingRast <- atan2(x=EstWestChng, y= NrthSthChng)*(180/pi)
# Turn according to the temporal trend  
  BearingRast1 <- app(c(BearingRast,
                        TempHetRast),
                      function(x){
                        ifelse(is.na(x[1]),
                               NA,
                               ifelse(c(x[2]<0),
                                      -x[1],x[1]))})    
  # make the bearing a value between 0 and 360
  BearingRast2 <- app(BearingRast1,
                      fun=function(x){ifelse(x<0,
                                             360+x,
                                             x)})
  # Turn these so that 0 bearing is Poleward
  BearingRast2a <- BearingRast2-90
  BearingRast2a <- app(BearingRast2a,
                       fun=function(x){ifelse(x<=0,
                                              360+x,
                                              x)})
  # Make the Bearing of areas that do not change 361DD
  BearingRast2a[][which(EstWestChng[]==0 & NrthSthChng[]==0)]<-361
  # Return the Bearing
  return(BearingRast2a)
}
#-------------------------------------------------------------------------------
## VelocityFnc: Function to estimate the velocity of a surface as the ratio 
##              between the spatial gradients and temporal changes.
VelocityFnc  <- function(TempHet,SpatHet){
  # Estimate the velocity of change
  Velocity <-  abs(TempHet)/SpatHet
  # Make the velocity for very slow areas (~1e-3km/year) zero
  Velocity[which(Velocity[]<0.001)] <- 0.001
  # Make the velocity for places with no temporal change zero  
  Velocity[which(TempHet[]==0)] <- 0.001
  # Make the velocity for flat areas 200km/year
  Velocity[which(SpatHet[]==0)] <- 200
  # Make the velocity for very fast areas (>200km/year) the max value 200km/year
  Velocity[which(Velocity[]>200)] <- 200
  return(Velocity)
}
#-------------------------------------------------------------------------------#-------------------------------------------------------------------------------
