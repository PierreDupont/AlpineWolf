## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------ IMPORT REQUIRED LIBRARIES ------
library(rgdal)
library(raster)
library(coda)
library(nimble)
library(nimbleSCR)
library(stringr)
library(abind)
library(R.utils)
library(sf)
library(fasterize)
library(dplyr)
library(lubridate)
library(stars)
library(RANN)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(gridExtra)


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.5.0"
thisDir <- file.path(analysisDir, modelName)
if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(!dir.exists(thisDir)){dir.create(thisDir)}
if(!dir.exists(file.path(thisDir, "input"))){dir.create(file.path(thisDir, "input"))}
if(!dir.exists(file.path(thisDir, "output"))){dir.create(file.path(thisDir, "output"))}

## LOAD HABITAT SPECIFICATIONS
load(file.path(analysisDir, modelName, "Habitat.RData"))

## LOAD DETECTORS SPECIFICATIONS
load(file.path(analysisDir, modelName, "Detectors.RData"))

## LOAD PROCESSED MCMC SAMPLES
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))

## LOAD MODEL INPUT 
load(file.path(thisDir, "input", paste0(modelName, "_1.RData")))

## LOAD POYGON OF ITALY AND NEIGHBOURING COUNTRIES
countries <- read_sf(file.path(dataDir,"GISData/Italy_borders/Italy_andBorderCountries_splitFrance.shp"))

## LOAD POLYGONS OF ITALY's REGIONS 
regions <- read_sf(file.path(dataDir,"GISData/Output_layout/Alpine_Regions.shp"))
regions <- st_transform(x = regions, crs = st_crs(countries))
regions$ID <- as.numeric(as.factor(regions$DEN_UTS))
plot(regions)

## LOAD POLYGON OF THE ITALIAN ALPS
alps <- read_sf(file.path(dataDir,"GISData/Output_layout/Italian_Alps.shp"))
alps <- st_transform(x = alps, crs = st_crs(countries))
plot(alps)


## -----------------------------------------------------------------------------
## ------ 1. CALCULATE AC-BASED DENSITY ------
## ------   1.1. PREPARE INPUT FOR EXTRACTION ------
##---- Disaggregate habitat raster to a 1x1km resolution
habitat.r <- disaggregate(x = habitat$raster, fact = 5)
italia.r <- disaggregate(x = habitat$Italia, fact = 5)

##---- Create a matrix of raster cellsIDs (including cells outside the habitat)
habitat.id <- matrix( data = 1:ncell(habitat.r),
                      nrow = nrow(habitat.r),
                      ncol = ncol(habitat.r),
                      byrow = TRUE)

##---- Rescale sxy coords to original projection and reproject to the new raster
dimnames(res$sims.list$s) <- list(1:dim(res$sims.list$s)[1],
                                  1:dim(res$sims.list$s)[2],
                                  c("x","y"))

s.original <- scaleCoordsToHabitatGrid(
  coordsData = res$sims.list$s,
  coordsHabitatGridCenter = coordinates(habitat$raster),
  scaleToGrid = F)$coordsDataScaled

##---- ... and reproject to the new raster
s.rescaled <- scaleCoordsToHabitatGrid(
  coordsData = s.original,
  coordsHabitatGridCenter = coordinates(habitat.r),
  scaleToGrid = T)$coordsDataScaled

##---- Create a matrix of Italian regions
##---- (rows == regions ; columns == habitat raster cells)
regions.r <- fasterize(sf = st_as_sf(regions),
                       raster = habitat.r,
                       field = "ID",
                       background = 0)
regions.r[regions.r[]==0] <- NA
regions.r <- regions.r + italia.r - 1
regions.unique <- sort(na.omit(unique(regions.r[])))
regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions.r[] == x}))
regions.rgmx[is.na(regions.rgmx)] <- 0
row.names(regions.rgmx) <- sort(regions$DEN_UTS)

##---- Create a matrix of the Alpine region
##---- (rows == regions ; columns == habitat raster cells)
alps.r <- fasterize( sf = st_as_sf(alps),
                     raster = habitat.r,
                     background = 0)
alps.r[alps.r[]==0] <- NA
alps.r <- alps.r + italia.r - 1
alps.rgmx <- matrix(alps.r[] == 1, nrow = 1)
alps.rgmx[is.na(alps.rgmx)] <- 0
row.names(alps.rgmx) <- "Italian Alps"

##---- Thin MCMC samples for faster calculation 
##---- (60000 iterations in too much for the getDensity fucntion at 1x1km)
iter <- seq(1,dim(res$sims.list$z)[1],10)



## ------   1.2. EXTRACT TOTAL & REGION-SPECIFIC DENSITIES ------
##---- Calculate total and regional densities
WA_Regions <- GetDensity(
  sx = s.rescaled[iter, ,1],
  sy = s.rescaled[iter, ,2],
  z = res$sims.list$z[iter, ],
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = F,
  regionID = regions.rgmx)

##---- Calculate sex and status-specific densities
WA_status <- list()
for(s in 0:1){
  WA_status[[s+1]] <- list()
  for(ss in 1:3){
    thisStatus <- (res$sims.list$z[iter, ] == 1) &
      (res$sims.list$sex[iter, ] == s) &
      (res$sims.list$status[iter, ] == ss)
    
    WA_status[[s+1]][[ss]] <- GetDensity(
      sx = s.rescaled[iter, ,1],
      sy = s.rescaled[iter, ,2],
      z = thisStatus,
      IDmx = habitat.id,
      aliveStates = 1,
      returnPosteriorCells = F,
      regionID = regions.rgmx)
  }#ss
}#s


## ------   1.3. EXTRACT ALPINE REGION DENSITY ------
WA_Alps <- GetDensity(
  sx = s.rescaled[iter, ,1],
  sy = s.rescaled[iter, ,2],
  z = res$sims.list$z[iter, ],
  IDmx = habitat.id,
  aliveStates = 1,
  returnPosteriorCells = F,
  regionID = alps.rgmx)


## ------   1.4. EXTRACT PREDICTED DENSITY ------
predDensities <- sapply(iter,
                        function(x){
                          intens <- c(exp(res$sims.list$betaHab[x, ] %*% t(nimData$hab.covs)))
                          pred <- res$sims.list$N[x]*(intens/sum(intens)) 
                          return(pred)
                        })

meanPred.r <- habitat$raster
meanPred.r[meanPred.r[] > 0] <- rowMeans(predDensities)
meanPred.r[is.na(habitat$Italia[])] <- NA



## ------   1.5. SAVE DENSITIES -----
save(WA_regions, WA_status, WA_Alps, predDensities, 
     file = file.path(thisDir, paste0(modelName, "_densities.RData")))



## -----------------------------------------------------------------------------
## ------ 2. REPORT .pdf ------
pdf(file = file.path(thisDir, paste0(modelName,"_Results.pdf")),
    width = 20, height = 15)

## ------   2.1. DENSITY MAPS ------

## ------     2.1.1. DENSITY MAPS ------
##---- Set color scale
maxDens <- max(WA_Regions$MeanCell)
cuts <- seq(0, maxDens, length.out = 100)   
colFunc <- colorRampPalette(c("white","slateblue","navyblue", "midnightblue"))#,"yellow","orange","red","red"
col <- colFunc(100)


##---- Plot overall density raster
par(mfrow = c(1,1), mar = c(4,4,0,0))

meanDensity.R <- habitat.r
meanDensity.R[ ] <- WA_Regions$MeanCell
meanDensity.R[is.na(italia.r[])] <- NA

plot( habitat$polygon, border = "white")
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( st_geometry(countries), add = T, lwd = 1)
plot( st_geometry(st_intersection(regions,countries)), add = T, lwd = 2)

mtext( text = paste( "N = ", round(WA_Regions$summary["Total",1],1),
                     " [", round(WA_Regions$summary["Total",4],1), " ; ",
                     round(WA_Regions$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)


##---- Plot Alpine region density raster
par(mfrow = c(1,1), mar = c(4,4,0,0))

alpDensity.R <- habitat.r
alpDensity.R[] <- WA_Alps$MeanCell
alpDensity.R[is.na(alps.r[])] <- NA

plot( habitat$polygon, border = "white")
plot( alpDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( st_geometry(countries), add = T, lwd = 1)
plot( st_geometry(st_intersection(alps,countries)), add = T, lwd = 2)

mtext( text = paste( "N = ", round(WA_Alps$summary["Total",1],1),
                     " [", round(WA_Alps$summary["Total",4],1), " ; ",
                     round(WA_Alps$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)



## ------  REALIZED vs. PREDICTED DENSITY  
par(mfrow = c(1,2), mar = c(1,1,1,1))
plot( habitat$polygon, border = "white")
mtext( text = "Realized density", side = 3, line = -4, font = 2)
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( st_geometry(countries), add = T, lwd = 1)
plot( st_geometry(st_intersection(regions,countries)), add = T, lwd = 2)

plot( habitat$polygon, border = "white")
mtext( text = "Predicted density", side = 3, line = -4, font = 2)
plot( meanPred.r/25, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( st_geometry(countries), add = T, lwd = 1)
plot( st_geometry(st_intersection(regions,countries)), add = T, lwd = 2)

## Add density legend
plot( meanDensity.R, breaks = cuts, col = col,
      legend.width = 1, legend.only = TRUE,
      axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 5),
                        labels = round(seq(0,maxDens,length.out = 5), digits = 2), 
                        cex.axis = 0.6),
      legend.args = list( text = 'Relative density (% of pop.km2)', line = -2,
                          side = 4, font = 2, cex = 0.8))



# ## ------  DETECTION CENTROIDS
# par(mfrow = c(1,1), mar = c(4,4,4,4))
# 
# ##-- Set color scale
# maxDens <- max(meanDensity.R[], na.rm = T)
# cuts <- seq(0, maxDens, length.out = 100)   
# colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
# col <- colFunc(100)
# 
# centroids <- st_drop_geometry(ngs) %>%
#   group_by(Genotype.ID) %>%
#   summarise(X = mean(CoordX),               ## Get total length searched in each detector grid cell
#             Y =  mean(CoordY)) 
# coordinates(centroids) <- cbind.data.frame(centroids$X,
#                                            centroids$Y)
# proj4string(centroids) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
# centroids <- st_as_sf(centroids)
# 
# plot( habitat$polygon, border = grey(0.3))
# plot( meanDensity.R, add = T,
#       breaks = cuts, col = col,
#       axes = F, box = F, bty = "n", legend = F)
# plot(  habitat$polygon, add = T)
# plot(st_geometry(countries), add = T, lwd=2)
# plot(centroids,
#      add = T, pch = 21, cex = 1.2,
#      bg = adjustcolor("red",0.2),
#      col = "red")
# mtext( text = "Wolf density with\ndetected individuals' centroid",
#        side = 3, font = 2, cex = 2)



## ------     2.1.2. DENSITY EFFECT PLOTS ------
par(mfrow = c(2,2), mar = c(6,6,0,0))
covNames <- names(nimData$hab.covs)
pred.hab.covs <- apply(nimData$hab.covs,
                       2,
                       function(c){
                         seq(min(c),
                             max(c),
                             length.out = 100)})
mean.hab.covs <- apply(nimData$hab.covs, 2, mean)
cols <- hcl.colors(length(mean.hab.covs))

for(b in 1:ncol(res$sims.list$betaHab)){
  intensity <- do.call( rbind,
                        lapply(1:length(res$sims.list$betaHab[ ,b]),
                               function(x){
                                 exp(res$sims.list$betaHab[x,b] * pred.hab.covs[ ,b] + 
                                       sum(res$sims.list$betaHab[x,-b] * mean.hab.covs[-b]))
                               }))
  
  
  mean.int <- colMeans(intensity)
  quant.int <- apply(intensity, 2, function(x)quantile(x,c(0.025,0.5,0.975)))
  maxD <- round(max(quant.int),3)
  plot( x = pred.hab.covs[ ,b],
        y = quant.int[2, ],
        type = "n", ylim = c(0, maxD), xlim = range(pred.hab.covs[ ,b]),
        ylab = "Density", xlab = covNames[b], axes = FALSE)
  minCov <- min(st_drop_geometry(habitat$grid[ ,covNames[b]]))
  maxCov <- max(st_drop_geometry(habitat$grid[ ,covNames[b]]))
  xLabels <- round(seq(minCov, maxCov, length.out = 10),2)
  axis(1,
       at = round(seq(min(pred.hab.covs[ ,b]), max(pred.hab.covs[ ,b]), length.out = 10),3),
       labels = xLabels, cex = 2,
       tck = 0.01, las = 1, hadj = 0.5)
  axis(2, at = seq(0,maxD,length.out = 6),
       labels = seq(0,maxD,length.out = 6),
       tck = 0.01, las = 1, cex = 2)
  
  polygon(x = c(pred.hab.covs[ ,b],rev(pred.hab.covs[ ,b])),
          y = c(quant.int[1, ],rev(quant.int[3, ])), border = F,
          col = adjustcolor(cols[b], alpha.f = 0.5))
  
  points( x = pred.hab.covs[ ,b],
          y = quant.int[2,],
          lwd = 2, type = "l", col = cols[b])
  
}
graphics.off()




## ------   2.2. DETECTION ------
pdf(file = file.path(thisDir, paste0(modelName,"_Detection.pdf")),
    width = 15, height = 15)

## ------     2.2.1. DETECTION MAPS ------
sex <- c("female","male")
status <- c("alpha","pup","other")
for(ss in 1:2){
  for(s in 1:3){
    detectors$grid[ ,paste0("p0_",sex[ss],"_",status[s])] <- c(ilogit(logit(res$mean$p0[s,ss]) +
                                                                        res$mean$betaDet %*% t(nimData$det.covs)))
  }
}
plot(detectors$grid[ ,c("p0_female_alpha", 
                        "p0_male_alpha",
                        "p0_female_pup",
                        "p0_male_pup",
                        "p0_female_other",
                        "p0_male_other")],
     key.pos = 4, breaks = seq(0,1,0.05)) 



## ------     2.2.2. DETECTION EFFECT PLOTS ------
par(mfrow = c(3,2), mar = c(8,8,0,4))
names(nimData$det.covs)

pred.det.covs <- apply(nimData$det.covs,
                       2,
                       function(c){
                         seq(min(c),
                             max(c),
                             length.out = 100)})
mean.det.covs <- apply(nimData$det.covs, 2, mean)
cols <- hcl.colors(length(mean.det.covs))

for(b in 1:ncol(res$sims.list$betaDet)){
  p0 <- do.call( rbind,
                 lapply(1:length(res$sims.list$betaDet[,b]),
                        function(x){
                          ilogit(logit(res$sims.list$p0[x,1,1]) +
                                   res$sims.list$betaDet[x,b] * pred.det.covs[ ,b]+ 
                                   sum(res$sims.list$betaDet[x,-b] * mean.det.covs[-b]))
                        }))
  quant.int <- apply(p0, 2, function(x)quantile(x,c(0.025,0.5,0.975)))
  maxp0 <- round(max(quant.int),3)
  plot( x = pred.det.covs[ ,b],
        y = quant.int[2, ],
        type = "n", ylim = c(0, maxp0), xlim = range(pred.det.covs[ ,b]),
        ylab = "p0", xlab = covNames[b], axes = FALSE)
  minCov <- min(st_drop_geometry(detectors$grid[ ,covNames[b]]))
  maxCov <- max(st_drop_geometry(detectors$grid[ ,covNames[b]]))
  xLabels <- round(seq(minCov, maxCov, length.out = 10),2)
  axis(1,
       at = round(seq(min(pred.det.covs[ ,b]), max(pred.det.covs[ ,b]), length.out = 10),3),
       labels = xLabels, cex = 3,
       tck = 0.01, las = 1, hadj = 0.5)
  axis(2, at = seq(0,maxp0,length.out = 6),
       labels = seq(0,maxp0,length.out = 6),
       tck = 0.01, las = 1, cex = 3)
  
  polygon(x = c(pred.det.covs[ ,b],rev(pred.det.covs[ ,b])),
          y = c(quant.int[1, ],rev(quant.int[3, ])), border = F,
          col = adjustcolor(cols[b], alpha.f = 0.5))
  
  points( x = pred.det.covs[ ,b],
          y = quant.int[2,],
          lwd = 2, type = "l", col = cols[b])
  
}


## ------   2.3. PARAMETER TABLES ------
## ------     2.3.1. SCR PARAMETERS ------
# paramSimple <- sapply(strsplit(colnames(res$sims.list), split = '\\['), '[', 1)
params.simple <- names(res$mean)[!names(res$mean) %in% c("s","z")]

params.means <- do.call(c, lapply(params.simple, function(x)res$mean[[x]]))
params.sd <- do.call(c, lapply(params.simple, function(x)res$sd[[x]]))
params.q2.5 <- do.call(c, lapply(params.simple, function(x)res$q2.5[[x]]))
params.q97.5 <- do.call(c, lapply(params.simple, function(x)res$q97.5[[x]]))
params.Rhat <- do.call(c, lapply(params.simple, function(x)res$Rhat[[x]]))
params.n.eff <- do.call(c, lapply(params.simple, function(x)res$n.eff[[x]]))

params.summary <- cbind.data.frame( params.means, params.sd,
                                    params.q2.5, params.q97.5, params.Rhat,
                                    params.n.eff)

names(params.summary) <- c("mean", "sd", "2.5%CI", "97.5%CI", "Rhat", "n.eff")

## WRITE TABLE
# write.csv( x = params.summary,
#            file = file.path(myVars$WD, myVars$modelName, "TABLES",
#                             paste( myVars$modelName, "_params.csv", sep = "")))
# 
# print( xtable(params.summary, type = "latex"),
#        floating = FALSE,# scalebox=.8,
#        add.to.row = list(list(seq(1, nrow(params.summary), by = 2)),"\\rowcolor[gray]{.95} "),
#        file = file.path(myVars$WD, myVars$modelName,"TABLES",
#                         paste( myVars$modelName, "_params.tex", sep = "")))




## ------     2.3.2. ABUNDANCES ------
WA_regions$summary


## ------  REGIONAL ABUNDANCES
##-- Export as .pdf
par(mfrow = c(1,1))
plot.new()
grid.table(round(WA_regions$summary,1))
mtext(text = "Abundance estimates\nby region",
      side = 3, outer = T, line = -20, cex = 2, font = 2)

sex <- c("female","male")
status <- c("alpha","pups","others")
par(mfrow = c(2,3))
for(s in 1:length(sex)){
  for(ss in 1:length(status)){
    par(mfrow = c(1,1))
    plot.new()
    grid.table(round(WA_status[[s]][[ss]]$summary,1))
    mtext(text = paste0("Abundance estimates\nfor ", sex[s], "_", status[ss]),
          side = 3, outer = T, line = -10, font = 2)
  }
}


## ------   2.4. SPACE-USE ------
# ## Save input files for a few iterations
# iter <- seq(2001, dim(res$sims.list$s)[1], by = 100)
# 
# WA_SpaceUse <- GetSpaceUse(  sx = res$sims.list$s[20, ,1],
#                              sy = res$sims.list$s[20, ,2],
#                              z = res$sims.list$z[20, ],
#                              sigma = matrix(res$sims.list$sigma[20,1,1],
#                                             length(iter),
#                                             ncol(res$sims.list$z)),
#                              habitatxy = habitat$scaledCoords,
#                              aliveStates = 1,
#                              regionID = it.rgmx,
#                              returnPosteriorCells = F)
# WA_SpaceUse$summary
# 
# ##-- Set color scale
# maxDens <- max(c(WA_Density$MeanCell, WA_SpaceUse$MeanCell))
# cuts <- seq(0, maxDens, length.out = 100)   
# colFunc <- colorRampPalette(c("white","slateblue","yellow","orange","red","red"))
# col <- colFunc(100)
# 
# 
# ## ------       5.1. AC-based DENSITY MAP 
# ital.R <- habitat.r
# ital.R[ ] <- WA_Density$MeanCell
# ital.R[is.na(habitat$Italia[])] <- NA
# 
# plot( habitat$polygon, col = "gray80", border = grey(0.3))
# plot( ital.R, add = T,
#       breaks = cuts, col = col,
#       axes = F, box = F, bty = "n", legend = F)
# plot( habitat$polygon, add = T, border = grey(0.3))
# plot( st_geometry(countries), add = T, lwd = 2)
# 
# mtext( text = paste( "N = ", round(WA_regions$summary["Total",1],1),
#                      " [", round(WA_regions$summary["Total",4],1), " ; ",
#                      round(WA_regions$summary["Total",5],1), "]", sep = ""),
#        side = 1, font = 2, cex = 1.5)
# 
# plot( meanDensity.R, breaks = cuts, col = col,
#       legend.width = 1, legend.only = TRUE,
#       axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 1),
#                         labels = round(seq(0,maxDens,length.out = 5), digits = 1), 
#                         cex.axis = 0.6),
#       legend.args = list( text = 'Density (ids.km2)', line = -2,
#                           side = 4, font = 2, cex = 0.8))
# 
# 
# 
# 
# ## ------       5.2. UD-based DENSITY MAP 
# habitat.r <- Habitat.list$habitat.r
# habitat.r[habitat.r == 0] <- NA
# admin.r <- rasterize(Cleanpoly, habitat.r)
# 
# meanDensity.R <- habitat.r
# meanDensity.R[meanDensity.R[]==1] <- WA_SpaceUse$MeanCell
# 
# 
# plot( habitat$polygon, col = "gray80", border = grey(0.3))
# plot( meanDensity.R, add = T,
#       breaks = cuts, col = col,
#       axes = F, box = F, bty = "n", legend = F)
# plot( habitat$polygon, add = T, border = grey(0.3))
# plot( st_geometry(countries), add = T, lwd = 2)
# 
# mtext( text = paste( "N = ", round(WA_SpaceUse$summary["Total",1],1),
#                      " [", round(WA_SpaceUse$summary["Total",4],1), " ; ",
#                      round(WA_SpaceUse$summary["Total",5],1), "]", sep = ""),
#        side = 1, font = 2, cex = 1.5)
# par(mfrow = c(1,1), mar = c(4,4,0,0))
# plot( Cleanpoly, col = "gray80", border = grey(0.3))
# plot( meanDensity.R, add = T,
#       breaks = cuts, col = col, colNA = NA,
#       axes = F, box = F, bty = "n", legend = F)
# plot( Cleanpoly, add = T, border = grey(0.1), lwd = 3)
# mtext( text = paste( "N = ", round(SpaceUse$summary["Total",1],1),
#                      " [", round(SpaceUse$summary["Total",4],1), " ; ",
#                      round(SpaceUse$summary["Total",5],1), "]", sep = ""),
#        side = 1, font = 2, cex = 1.5)
# mtext( text = "Total utilization distribution", side = 3, font = 2, cex = 2, line = -2)
# 
# ##-- Add density legend
# plot( meanDensity.R, breaks = cuts, col = col,
#       legend.width = 1, legend.only = TRUE,
#       axis.args = list( at = round(seq(0,maxDens,length.out = 5), digits = 1),
#                         labels = round(seq(0,maxDens,length.out = 5), digits = 1), 
#                         cex.axis = 0.6),
#       legend.args = list( text = 'Density (ids.km2)', line = -2,
#                           side = 4, font = 2, cex = 0.8))

graphics.off()



## -----------------------------------------------------------------------------