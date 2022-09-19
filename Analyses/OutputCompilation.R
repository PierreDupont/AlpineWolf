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
library(MetBrewer)


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.5.2"
thisDir <- file.path(analysisDir, modelName)
if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(!dir.exists(thisDir)){dir.create(thisDir)}
if(!dir.exists(file.path(thisDir, "input"))){dir.create(file.path(thisDir, "input"))}
if(!dir.exists(file.path(thisDir, "output"))){dir.create(file.path(thisDir, "output"))}

## LOAD HABITAT SPECIFICATIONS
load(file.path(analysisDir, modelName, "Habitat.RData"))

## LOAD DETECTORS SPECIFICATIONS
load(file.path(analysisDir, modelName, "Detectors.RData"))

# LOAD PROCESSED MCMC SAMPLES
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))

## LOAD MODEL INPUT 
load(file.path(thisDir, "input", paste0(modelName, "_1.RData")))

## LOAD POLYGON OF ITALY AND NEIGHBOURING COUNTRIES
countries <- read_sf(file.path(dataDir,"GISData/Italy_borders/Italy_andBorderCountries_splitFrance.shp"))

## LOAD POLYGONS OF ITALY's REGIONS 
regions <- read_sf(file.path(dataDir,"GISData/Output_layout/Alpine_Regions.shp"))
regions <- st_transform(x = regions, crs = st_crs(countries))
regions$ID <- as.numeric(as.factor(regions$DEN_UTS))

## LOAD POLYGON OF THE ITALIAN ALPS
alps <- read_sf(file.path(dataDir,"GISData/Output_layout/Italian_Alps.shp"))
alps <- st_transform(x = alps, crs = st_crs(countries))

# ## LOAD POLYGON OF THE SEARCHED AREAS
# presTrans <- read_sf(file.path(dataDir,"GISData/Merged_PRESENCE_TRANSECTS/Presence_transects10km_buffer.shp"))
# presTrans <- st_intersection(presTrans, regions)
# presTrans <- st_cast(presTrans,"MULTIPOLYGON")

## LOAD POLYGON OF THE SEARCHED AREAS
east <- read_sf(file.path(dataDir,"GISData/Output_layout/East_tot.shp"))
west <- read_sf(file.path(dataDir,"GISData/Output_layout/West_tot.shp"))

## LOAD POLYGON OF THE STUDY AREA
studyAreaGrid <- read_sf(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- st_transform(x = studyAreaGrid, crs = st_crs(countries))
studyArea <- studyAreaGrid %>%
  st_snap(x = ., y = ., tolerance = 0.0001) %>%
  st_union() 

## -----------------------------------------------------------------------------
## ------ 1. CALCULATE AC-BASED DENSITY ------
# ## ------   1.1. PREPARE INPUT FOR EXTRACTION ------
# ##---- Disaggregate habitat raster to a 1x1km resolution
# habitat.r <- disaggregate(x = habitat$raster, fact = 1)
# italia.r <- disaggregate(x = habitat$Italia, fact = 1)
# 
# ##---- Create a matrix of raster cellsIDs (including cells outside the habitat)
# habitat.id <- matrix( data = 1:ncell(habitat.r),
#                       nrow = nrow(habitat.r),
#                       ncol = ncol(habitat.r),
#                       byrow = TRUE)
# 
# ##---- Rescale sxy coords to original projection and reproject to the new raster
# dimnames(res$sims.list$s) <- list(1:dim(res$sims.list$s)[1],
#                                   1:dim(res$sims.list$s)[2],
#                                   c("x","y"))
# s.original <- scaleCoordsToHabitatGrid(
#   coordsData = res$sims.list$s,
#   coordsHabitatGridCenter = coordinates(habitat$raster),
#   scaleToGrid = F)$coordsDataScaled
# 
# ##---- ... and reproject to the new raster
# s.rescaled <- scaleCoordsToHabitatGrid(
#   coordsData = s.original,
#   coordsHabitatGridCenter = coordinates(habitat.r),
#   scaleToGrid = T)$coordsDataScaled
# 
# ##---- Create a matrix of East/west regions
# ##---- (rows == regions ; columns == habitat raster cells)
# regions.r <- fasterize(sf = st_as_sf(regions),
#                        raster = habitat.r,
#                        field = "ID",
#                        background = NA)
# west <- st_buffer(st_union(west), dist = 10000)
# west <- st_difference(west,st_union(east))
# west.r <- fasterize( sf = st_as_sf(west),
#                      raster = habitat.r,
#                      background = 0)
# reg.r <- regions.r
# reg.r[reg.r>0] <- 1
# regions2.r <- reg.r + west.r 
# regions2.r[regions2.r[]==1] <- "east"
# regions2.r[regions2.r[]==2] <- "west"
# regions.r <- regions.r + italia.r - 1
# 
# regions.unique <- sort(na.omit(unique(regions2.r[])))
# regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions2.r[] == x}))
# regions.rgmx[is.na(regions.rgmx)] <- 0
# row.names(regions.rgmx) <- regions.unique
# 
# 
# # ##---- Create a matrix of Italian regions
# # ##---- (rows == regions ; columns == habitat raster cells)
# # regions.r <- fasterize(sf = presTrans,
# #                        raster = habitat.r,
# #                        field = "ID",
# #                        background = 0)
# # regions.r[regions.r[]==0] <- NA
# # regions.r <- regions.r + italia.r - 1
# # for(i in 1:length(presTrans$DEN_UTS)){
# #   regions.r[regions.r %in% regions$ID[i]] <- presTrans$DEN_UTS[i]
# # }
# # regions.unique <- sort(na.omit(unique(regions.r[])))
# # regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions.r[] == x}))
# # regions.rgmx[is.na(regions.rgmx)] <- 0
# # row.names(regions.rgmx) <- regions.unique
# 
# 
# ##---- Create a matrix of the Alpine region
# ##---- (rows == regions ; columns == habitat raster cells)
# alps.r <- fasterize( sf = st_as_sf(alps),
#                      raster = habitat.r,
#                      background = NA)
# alps.r[alps.r[]==0] <- NA
# alps.r <- alps.r + italia.r - 1
# alps.rgmx <- matrix(alps.r[] == 1, nrow = 1)
# alps.rgmx[is.na(alps.rgmx)] <- 0
# row.names(alps.rgmx) <- "Italian Alps"
# 
# ##---- Thin MCMC samples for faster calculation 
# ##---- (60000 iterations in too much for the getDensity fucntion at 1x1km)
# iter <- seq(1,dim(res$sims.list$z)[1],length.out = 20000)
# 
# 
# 
# ## ------   1.2. EXTRACT TOTAL & REGION-SPECIFIC DENSITIES ------
# ##---- Calculate total and regional densities
# WA_regions <- GetDensity(
#   sx = s.rescaled[iter, ,1],
#   sy = s.rescaled[iter, ,2],
#   z = res$sims.list$z[iter, ],
#   IDmx = habitat.id,
#   aliveStates = 1,
#   returnPosteriorCells = F,
#   regionID = regions.rgmx)
# 
# ##---- Calculate sex and status-specific densities
# WA_status <- list()
# for(s in 0:1){
#   WA_status[[s+1]] <- list()
#   for(ss in 1:3){
#     thisStatus <- (res$sims.list$z[iter, ] == 1) &
#       (res$sims.list$sex[iter, ] == s) &
#       (res$sims.list$status[iter, ] == ss) 
#     
#     WA_status[[s+1]][[ss]] <- GetDensity(
#       sx = s.rescaled[iter, ,1],
#       sy = s.rescaled[iter, ,2],
#       z = thisStatus,
#       IDmx = habitat.id,
#       aliveStates = 1,
#       returnPosteriorCells = F,
#       regionID = regions.rgmx)
#   }#ss
# }#s
# 
# ##---- Name the object to fill up rapidly the tables 
# names(WA_status) <- c("F","M")
# names(WA_status[["F"]]) <- c("Alpha", "Pups", "Others")
# names(WA_status[["M"]]) <- c("Alpha", "Pups", "Others")
# 
# 
# ##---- Calculate sex and status-specific densities
# augmented <- (sum(nimData$y[ ,1] > 0)+1):nimConstants$n.individuals
# 
# WA_undetected <- list()
# for(s in 0:1){
#   WA_undetected[[s+1]] <- list()
#   for(ss in 1:3){
#     thisStatus <- (res$sims.list$z[iter,augmented] == 1) &
#       (res$sims.list$sex[iter,augmented] == s) &
#       (res$sims.list$status[iter,augmented] == ss) 
#     
#       WA_undetected[[s+1]][[ss]] <- GetDensity(
#       sx = s.rescaled[iter,augmented,1],
#       sy = s.rescaled[iter,augmented,2],
#       z = thisStatus,
#       IDmx = habitat.id,
#       aliveStates = 1,
#       returnPosteriorCells = F,
#       regionID = regions.rgmx)
#   }#ss
# }#s
# 
# ##---- Name the object to fill up rapidly the tables 
# names(WA_undetected) <- c("F","M")
# names(WA_undetected[["F"]]) <- c("Alpha", "Pups", "Others")
# names(WA_undetected[["M"]]) <- c("Alpha", "Pups", "Others")
# 
# 
# 
# ## ------   1.3. EXTRACT ALPINE REGION DENSITY ------
# WA_alps <- GetDensity(
#   sx = s.rescaled[iter, ,1],
#   sy = s.rescaled[iter, ,2],
#   z = res$sims.list$z[iter, ],
#   IDmx = habitat.id,
#   aliveStates = 1,
#   returnPosteriorCells = F,
#   regionID = alps.rgmx)
# 
# 
# ## ------   1.4. EXTRACT PREDICTED DENSITY ------
# predDensities <- sapply(iter,
#                         function(x){
#                           intens <- c(exp(res$sims.list$betaHab[x, ] %*% t(nimData$hab.covs)))
#                           pred <- res$sims.list$N[x]*(intens/sum(intens))
#                           return(pred)
#                         })
# 
# meanPred.r <- habitat$raster
# meanPred.r[meanPred.r[] > 0] <- rowMeans(predDensities)
# meanPred.r[is.na(habitat$Italia[])] <- NA
# 
# 
# 
# ## ------   1.5. EXTRACT TOTAL DENSITY (for plots) ------
# ##---- Create a density map for plotting at a a 1x1km resolution
# habitat.r <- disaggregate(x = habitat$raster, fact = 5)
# italia.r <- disaggregate(x = habitat$Italia, fact = 5)
# s.original <- scaleCoordsToHabitatGrid(
#   coordsData = res$sims.list$s,
#   coordsHabitatGridCenter = coordinates(habitat$raster),
#   scaleToGrid = F)$coordsDataScaled
# s.rescaled <- scaleCoordsToHabitatGrid(
#   coordsData = s.original,
#   coordsHabitatGridCenter = coordinates(habitat.r),
#   scaleToGrid = T)$coordsDataScaled
# regions.r <- fasterize( sf = st_as_sf(regions),
#                         raster = habitat.r,
#                         background = NA)
# regions.r[regions.r[]==0] <- NA
# regions.r <- regions.r + italia.r - 1
# regions.rgmx <- matrix(regions.r[] == 1, nrow = 1)
# regions.rgmx[is.na(regions.rgmx)] <- 0
# row.names(regions.rgmx) <- "Italian Alps"
# habitat.id <- matrix( data = 1:ncell(habitat.r),
#                       nrow = nrow(habitat.r),
#                       ncol = ncol(habitat.r),
#                       byrow = TRUE)
# iter <- seq(1,dim(res$sims.list$z)[1],length.out = 10000)
# WA_plot <- GetDensity(
#   sx = s.rescaled[iter, ,1],
#   sy = s.rescaled[iter, ,2],
#   z = res$sims.list$z[iter, ],
#   IDmx = habitat.id,
#   aliveStates = 1,
#   returnPosteriorCells = F,
#   regionID = regions.rgmx)
# 
# 
# 
# ## ------   1.6. SAVE DENSITIES -----
# save(WA_regions, WA_status, WA_alps, meanPred.r,  WA_undetected,
#      file = file.path(thisDir, paste0(modelName, "_densities.RData")))

load(file = file.path(thisDir, paste0(modelName, "_densities.RData")))


## -----------------------------------------------------------------------------
## ------ 2. REPORT .pdf ------
pdf(file = file.path(thisDir, paste0(modelName,"_Results.pdf")),
    width = 20, height = 15)

## ------   2.1. DENSITY MAPS ------
## ------     2.1.1. DENSITY MAPS ------
##---- Set color scale
maxDens <- max(WA_plot$MeanCell)
cuts <- seq(0, maxDens, length.out = 100)   
#cuts <- c(0,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256)   

colFunc <- colorRampPalette(c("black", "orange", "yellow", "white"))#,"yellow","orange","red","red"
col <- colFunc(100)
#col <- hcl.colors(10)

##---- Plot overall density raster
# pdf(file = file.path(thisDir, paste0(modelName,"_DensityMap.pdf")),
#     width = 20, height = 15)
par(mfrow = c(1,1), mar = c(4,4,0,0),bg="white")

habitat.r <- disaggregate(x = habitat$raster, fact = 5)
italia.r <- disaggregate(x = habitat$Italia, fact = 5)
meanDensity.R <- habitat.r
meanDensity.R[ ] <- WA_plot$MeanCell
meanDensity.R[is.na(italia.r[])] <- NA

bground <- st_buffer(countries[countries$CNTR_CODE == "IT", ], 50000)
bground <- st_difference(st_union(bground), st_union(countries))
plot(st_geometry(st_union(st_intersection(countries,studyArea))),
     col = "black", border = "black")
#plot( st_geometry(countries), col="black", border = "white", add = T, lwd = 1)
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( st_geometry(st_difference(countries,studyArea)),
      add = T, lwd = 2, col = "black",border = "white")
plot( st_geometry(bground),
      add = T, lwd = 2, col = "white",border = "white")
# plot( st_geometry(st_union(st_intersection(countries,studyArea))),
#       add = T, lwd = 2, border = "black",border="white")
# mtext( text = paste( "N = ", round(WA_regions$summary["Total",1],1),
#                      " [", round(WA_regions$summary["Total",4],1), " - ",
#                      round(WA_regions$summary["Total",5],1), "]", sep = ""),
#        side = 1, line = -5,font = 2, cex = 1.5, col = "white")
# dev.off()

##---- Plot Alpine region density raster
par(mfrow = c(1,1), mar = c(4,4,0,0))

alpDensity.R <- habitat$raster
alpDensity.R[] <- WA_alps$MeanCell
alpDensity.R[is.na(alps.r[])] <- NA

plot( habitat$polygon, border = "white")
plot( st_geometry(countries), col="black", border = "white", add = T, lwd = 1)
plot( alpDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( st_geometry(countries), add = T, lwd = 1, border = "white")
plot( st_geometry(st_intersection(alps,countries)), add = T, lwd = 1, border = "white")

mtext( text = paste( "N = ", round(WA_alps$summary["Total",1],1),
                     " [", round(WA_alps$summary["Total",4],1), " ; ",
                     round(WA_alps$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)


## ------  REALIZED vs. PREDICTED DENSITY  
par(mfrow = c(1,2), mar = c(1,1,1,1))
plot( habitat$polygon, border = "white")
plot( st_geometry(countries), col="black", border = "white", add = T, lwd = 1)
mtext( text = "Realized density", side = 3, line = -4, font = 2)
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( st_geometry(countries), add = T, lwd = 1, border = "white")

plot( habitat$polygon, border = "white")
mtext( text = "Predicted density", side = 3, line = -4, font = 2)
plot( st_geometry(countries), col="black", border = "white", add = T, lwd = 1)
plot( meanPred.r/25, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( st_geometry(countries), add = T, lwd = 1, border = "white")

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
par(mfrow = c(1,1), mar = c(10,10,2,8))
covNames <- names(nimData$hab.covs)
cols <- met.brewer(name="Isfahan1",n=length(covNames),type="discrete")
#cols <- hcl.colors(length(mean.hab.covs))

pred.hab.covs <- apply(nimData$hab.covs,
                       2,
                       function(c){
                         seq(min(c),
                             max(c),
                             length.out = 100)})
mean.hab.covs <- apply(nimData$hab.covs, 2, mean)

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
        y = quant.int[2, ], cex.lab = 3,
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




## ------   2.2. DETECTION ------
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
     key.pos = 4,
     breaks = c(0,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1)) 



## ------     2.2.2. DETECTION EFFECT PLOTS ------
# pdf(file = file.path(thisDir, paste0(modelName,"_DensityEffectPlots.pdf")),
#     width = 20, height = 15)
par(mfrow = c(1,1), mar = c(10,10,2,6))
covNames <- names(nimData$det.covs)
cols <- met.brewer(name="Isfahan1",n=length(covNames),type="discrete")

pred.det.covs <- apply(nimData$det.covs,
                       2,
                       function(c){
                         seq(min(c),
                             max(c),
                             length.out = 100)})
mean.det.covs <- apply(nimData$det.covs, 2, mean)

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
        ylab = "p0", xlab = covNames[b], axes = FALSE, cex.lab = 3)
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


## ------   2.3. SIGMA ------
par(mfrow = c(1,1), mar = c(15,25,15,25))
myCols <- matrix(met.brewer(name="Isfahan1",n=8,type="discrete")[c(2:4,7:5)],
                 nrow = 3, ncol=2, byrow = F)
valNoise <- c(-0.2,0,0.2)
plot(1,1,
     xlim = c(0.5,2.5),
     ylim = c(0,6.5),
     type = "n", xaxt = "n", main = "",
     ylab = expression(paste(sigma, " (km)")),
     xlab = "", axes = F)
axis(1, at = 1:2, labels = sex)
axis(2, at = seq(0,6.5,length.out = 6), labels = seq(0,6.5,length.out = 6)*5)
for(s in 1:length(sex)){
  for(ss in 1:length(status)){
    temp <- res$sims.list$sigma[ ,ss,s]
    
    plot.violins2( dat.list = list(temp),
                   x = s + valNoise[ss],
                   at = s + valNoise[ss],
                   add = T,
                   col = myCols[ss,s],
                   alpha = 0.9,
                   violin.width = 0.02,
                   border.col = myCols[ss,s],
                   scale.width = T)
  }#status
}#sex
legend(title = "Females",
       legend = c("Alpha","Pup","Other"),
       x = 0.5, y = 6.5, fill = myCols[1:3])
legend(title = "Males",
       legend = c("Alpha","Pup","Other"),
       x = 0.9, y = 6.5, fill = myCols[4:6])




## ------   2.4. PARAMETER TABLES ------
## ------     2.4.1. SCR PARAMETERS ------
sex <- c("Females","Males")
status <- c("Alpha","Pups","Others")

##---- Create empty table
popCompo  <-  matrix(NA, nrow = 2, ncol = 3)
colnames(popCompo) <- status
rownames(popCompo) <- sex

for(ss in 1:length(status)){
  for(s in 1:length(sex)){
    popCompo[s,ss] <- paste0(
      format(round( res$mean$theta[ss,s], 2), nsmall = 1),
      " (",  round(res$q2.5$theta[ss,s],2),
      "-", round(res$q97.5$theta[ss,s],2), ")")
  }#s
}#ss

##----- Export as .pdf
par(mfrow = c(1,1))
plot.new()
grid.table(popCompo)
mtext( text = "Alpine wolf population composition", 
       side = 3, outer = T, line = -10, font = 2)

# ##----- Export as .csv
# write.csv( abundanceTable,
#           file = file.path(thisDir, paste0(modelName, "_TableAbundance.csv")))
# 


## ------     2.4.2. ABUNDANCES PER REGION/SEX/STATUS ------
sex <- c("F","M")
status <- c("Alpha","Pups","Others")

##---- Create empty table
abundanceTable  <-  matrix(NA, nrow = 4, ncol = 12)
colnames(abundanceTable) <- c(rep("Alpha",3), rep("Pups",3), rep("Others",3), rep("Total",3))
rownames(abundanceTable) <- c("Region", row.names(WA_regions$summary))
abundanceTable[1,] <- c(rep(c("F", "M", "Total"), 4))

for(ss in 1:length(status)){
  colPos <- which(colnames(abundanceTable) %in% status[ss])
  for(s in 1:length(sex)){
    abundanceTable[2:nrow(abundanceTable),colPos[s]] <- paste0(
      format(round( WA_status[[sex[s]]][[status[ss]]]$summary[,"mean"],1), nsmall = 1),
      " (", WA_status[[sex[s]]][[status[ss]]]$summary[,"95%CILow"],
      "-", WA_status[[sex[s]]][[status[ss]]]$summary[,"95%CIHigh"], ")")
  }#s
  
  ##---- Calculate status-specific TOTALS from posteriors
  TotalTmp <- WA_status[[sex[1]]][[status[ss]]]$PosteriorRegions + WA_status[[sex[2]]][[status[ss]]]$PosteriorRegions
  abundanceTable[2:nrow(abundanceTable),colPos[3]] <- c(
    paste0( format(round(apply(TotalTmp, 1, mean), 1), nsmall = 1), " (",
            apply(TotalTmp,1,function(x) quantile(x,prob=c(0.025))), "-",
            apply(TotalTmp,1,function(x) quantile(x,prob=c(0.975))), ")"),
    ## TOTAL    
    paste0(format(round(mean(colSums(TotalTmp)), 1), nsmall = 1), " (",
           format(round(quantile(colSums(TotalTmp),prob=c(0.025)),0), nsmall = 0), "-",
           format(round(quantile(colSums(TotalTmp),prob=c(0.975)),0), nsmall = 0), ")"))
}#ss

##---- Get TOTALS per sex
colPos <- which(colnames(abundanceTable) %in% "Total")
for(s in 1:length(sex)){
  TotalTmp <-  WA_status[[sex[s]]][[1]]$PosteriorRegions + WA_status[[sex[s]]][[2]]$PosteriorRegions +
    WA_status[[sex[s]]][[3]]$PosteriorRegions
  abundanceTable[2:nrow(abundanceTable), colPos[s]] <- c(
    paste0(format(round(apply(TotalTmp,1,mean),1), nsmall = 1), " (",
           apply(TotalTmp,1,function(x) quantile(x,prob=c(0.025))), "-",
           apply(TotalTmp,1,function(x) quantile(x,prob=c(0.975))), ")"),
    ## TOTAL    
    paste0(format(round(mean(colSums(TotalTmp)),0), nsmall = 0), " (",
           quantile(colSums(TotalTmp),prob=c(0.025)), "-",
           quantile(colSums(TotalTmp),prob=c(0.975)), ")"))
}#s

##----- Get TOTAL per region both sex
abundanceTable[2:nrow(abundanceTable),colPos[3]] <- paste0(
  format(round(WA_regions$summary[,"mean"],1), nsmall = 1), " (",
  WA_regions$summary[,"95%CILow"], "-",
  WA_regions$summary[,"95%CIHigh"], ")")

##----- Export as .pdf
par(mfrow = c(1,1))
plot.new()
grid.table(abundanceTable)
mtext( text = "Alpine wolf abundance estimates", 
       side = 3, outer = T, line = -10, font = 2)

##----- Export as .csv
write.csv( abundanceTable,
           file = file.path(thisDir, paste0(modelName, "_TableAbundance.csv")))


dev.off()
## -----------------------------------------------------------------------------
## ------ 3. % of detected individuals ------
percentUndetected <- list()
for(s in 0:1){
  percentUndetected[[s+1]] <- list()
  for(ss in 1:3){
    percentUndetected[[s+1]][[ss]] <- WA_undetected[[s+1]][[ss]]$PosteriorAllRegions/WA_status[[s+1]][[ss]]$PosteriorAllRegions
  }#ss
  percentUndetected[[s+1]][[4]] <- (WA_undetected[[s+1]][[1]]$PosteriorAllRegions + 
                                    WA_undetected[[s+1]][[2]]$PosteriorAllRegions +
                                    WA_undetected[[s+1]][[3]]$PosteriorAllRegions)/ (WA_status[[s+1]][[1]]$PosteriorAllRegions +
                                                                                       WA_status[[s+1]][[2]]$PosteriorAllRegions +
                                                                                       WA_status[[s+1]][[3]]$PosteriorAllRegions)
}#s

percentUndetected[[3]] <- list()
for(ss in 1:3){
  percentUndetected[[3]][[ss]] <- (WA_undetected[[1]][[ss]]$PosteriorAllRegions + 
                                   WA_undetected[[2]][[ss]]$PosteriorAllRegions) / (WA_status[[1]][[ss]]$PosteriorAllRegions +
                                                                                      WA_status[[2]][[ss]]$PosteriorAllRegions)
}#ss
percentUndetected[[3]][[4]] <-  (WA_undetected[[1]][[1]]$PosteriorAllRegions + 
                                 WA_undetected[[1]][[2]]$PosteriorAllRegions +
                                 WA_undetected[[1]][[3]]$PosteriorAllRegions + 
                                 WA_undetected[[2]][[1]]$PosteriorAllRegions + 
                                 WA_undetected[[2]][[2]]$PosteriorAllRegions +
                                 WA_undetected[[2]][[3]]$PosteriorAllRegions) / (WA_status[[1]][[1]]$PosteriorAllRegions +
                                                                                    WA_status[[1]][[2]]$PosteriorAllRegions +
                                                                                    WA_status[[1]][[3]]$PosteriorAllRegions +
                                                                                    WA_status[[2]][[1]]$PosteriorAllRegions +
                                                                                    WA_status[[2]][[2]]$PosteriorAllRegions +
                                                                                    WA_status[[2]][[3]]$PosteriorAllRegions)

##---- Name the object to fill up rapidly the tables
names(percentUndetected) <- c("F","M","Total")
names(percentUndetected[["F"]]) <- c("Alpha", "Pups", "Others", "Total")
names(percentUndetected[["M"]]) <- c("Alpha", "Pups", "Others", "Total")
names(percentUndetected[["Total"]]) <- c("Alpha", "Pups", "Others", "Total")


mean(percentUndetected$F$Alpha)
mean(percentUndetected$F$Pups)
mean(percentUndetected$F$Others)
mean(percentUndetected$F$Total)

mean(percentUndetected$M$Alpha)
mean(percentUndetected$M$Pups)
mean(percentUndetected$M$Others)
mean(percentUndetected$M$Total)

mean(percentUndetected$Total$Alpha)
mean(percentUndetected$Total$Pups)
mean(percentUndetected$Total$Others)
mean(percentUndetected$Total$Total)

quantile(percentUndetected$Total$Alpha, probs = c(0.025,0.975))
quantile(percentUndetected$Total$Pups, probs = c(0.025,0.975))
quantile(percentUndetected$Total$Others, probs = c(0.025,0.975))
quantile(percentUndetected$Total$Total, probs = c(0.025,0.975))


