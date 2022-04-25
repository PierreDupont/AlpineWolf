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
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.5.2"
thisDir <- file.path(analysisDir, modelName)

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

## LOAD POLYGON OF THE SEARCHED AREAS
presTrans <- read_sf(file.path(dataDir,"GISData/Merged_PRESENCE_TRANSECTS/Presence_transects10km_buffer.shp"))
presTrans <- st_intersection(presTrans, regions)
presTrans <- st_cast(presTrans,"MULTIPOLYGON")

east <- read_sf(file.path(dataDir,"GISData/Output_layout/East_tot.shp"))
west <- read_sf(file.path(dataDir,"GISData/Output_layout/West_tot.shp"))

##--- Study area grid
studyAreaGrid <- read_sf(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- st_transform(x = studyAreaGrid, crs = st_crs(countries))
studyArea <- studyAreaGrid %>%
  st_snap(x = ., y = ., tolerance = 0.0001) %>%
  st_union() %>%
  st_intersection(.,countries)

##---- Load GPS search transects
transects <- read_sf(file.path(dataDir,"GISData/Transects_Wolfalps20202021/paths_completeness/paths_completeness.shp"))




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
                       background = NA)
# # regions.r[regions.r[]==0] <- NA
# regions.r <- regions.r + italia.r - 1
# for(i in 1:length(regions$DEN_UTS)){
#   regions.r[regions.r %in% regions$ID[i]] <- regions$DEN_UTS[i]
# }
# regions.unique <- sort(na.omit(unique(regions.r[])))
# regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions.r[] == x}))
# regions.rgmx[is.na(regions.rgmx)] <- 0
# row.names(regions.rgmx) <- regions.unique


##---- Create a matrix of East/West regions
##---- (rows == regions ; columns == habitat raster cells)
west <- read_sf(file.path(dataDir,"GISData/Output_layout/West_tot.shp"))
west <- st_buffer(st_union(west), dist = 10000)
west <- st_difference(west,st_union(east))
west.r <- fasterize( sf = st_as_sf(west),
                     raster = habitat.r,
                     background = 0)
reg.r <- regions.r
reg.r[reg.r>0] <- 1
regions2.r <- reg.r + west.r 
regions2.r[regions2.r[]==1] <- "east"
regions2.r[regions2.r[]==2] <- "west"

regions.unique <- sort(na.omit(unique(regions2.r[])))
regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions2.r[] == x}))
regions.rgmx[is.na(regions.rgmx)] <- 0
row.names(regions.rgmx) <- regions.unique



##---- Create a matrix of Italian regions
##---- (rows == regions ; columns == habitat raster cells)
regions.r <- fasterize(sf = st_as_sf(regions),
                       raster = habitat.r,
                       field = "ID",
                       background = NA)
# regions.r[regions.r[]==0] <- NA
regions.r <- regions.r + italia.r - 1
# for(i in 1:length(regions$DEN_UTS)){
#   regions.r[regions.r %in% regions$ID[i]] <- regions$DEN_UTS[i]
# }
regions.unique <- sort(na.omit(unique(regions.r[])))
regions.rgmx <- do.call(rbind, lapply(regions.unique, function(x){regions.r[] == x}))
regions.rgmx[is.na(regions.rgmx)] <- 0
row.names(regions.rgmx) <- regions.unique

##---- Create a matrix of the Alpine region
##---- (rows == regions ; columns == habitat raster cells)
alps.r <- fasterize( sf = st_as_sf(alps),
                     raster = habitat.r,
                     background = NA)
alps.r[alps.r[]==0] <- NA
alps.r <- alps.r + italia.r - 1
alps.rgmx <- matrix(alps.r[] == 1, nrow = 1)
alps.rgmx[is.na(alps.rgmx)] <- 0
row.names(alps.rgmx) <- "Italian Alps"

##---- Thin MCMC samples for faster calculation 
##---- (60000 iterations in too much for the getDensity fucntion at 1x1km)
iter <- seq(1,dim(res$sims.list$z)[1],length.out = 2000)



## ------   1.2. EXTRACT TOTAL & REGION-SPECIFIC DENSITIES ------
##---- Calculate total and regional densities
WA_regions <- GetDensity(
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

##---- Name the object to fill up rapidly the tables 
names(WA_status)<- c("F","M")
names(WA_status[["F"]]) <- c("Alpha", "Pups", "Others")
names(WA_status[["M"]]) <- c("Alpha", "Pups", "Others")



## ------   1.3. EXTRACT ALPINE REGION DENSITY ------
WA_alps <- GetDensity(
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


# ## ------   1.5. SAVE DENSITIES -----
# save(WA_regions, WA_status, WA_alps, meanPred.r,
#      file = file.path(thisDir, paste0(modelName, "_densities.RData")))

load(file = file.path(thisDir, paste0(modelName, "_densities.RData")))



## -----------------------------------------------------------------------------
## ------ 2. REPORT FIGURES ------
## ------   2.1. DENSITY MAPS ------
## ------     2.1.1. DENSITY MAPS ------
##---- Plot overall density raster
par(mfrow = c(1,1), mar = c(4,4,0,0), bg = "white")
meanDensity.R <- habitat.r
meanDensity.R[ ] <- WA_regions$MeanCell
meanDensity.R[is.na(italia.r[])] <- NA
plot(meanDensity.R)

##---- Smooth using a moving window 
f.rast <- function(x) ifelse(is.na(x[13]), NA, mean(x,na.rm = T)) 
MovingWindowSize <-  matrix(1,5,5)
meanDensity.R <- focal(meanDensity.R, MovingWindowSize, f.rast)
plot(meanDensity.R)



##---- Plot 1
##---- Set color scale
maxDens <- max(meanDensity.R[],na.rm = T)#max(WA_regions$MeanCell)
cuts <- seq(0, maxDens, length.out = 100)   


colFunc <- colorRampPalette(c("gray80", rev(terrain.colors(12))))
col <- colFunc(100)
{
pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_map1.pdf" )),
    width = 10, height = 8)
par(mfrow = c(1,1), mar = c(4,4,0,0),bg="white")
plot(st_geometry(st_union(st_intersection(countries,studyArea))),
     col = "white", border = "black")
plot( st_geometry(countries), add = T, lwd = 1, border = "black", col = "gray80")
plot( meanDensity.R, add = T,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot(st_geometry(st_union(st_intersection(countries,studyArea))),
     add=T, border = "black")
graphics.off()


pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_legend1.pdf" )),
    width = 1.75, height = 5)
par(mar=c(1,1,1,4),las=1,cex.axis=1.6,xaxs="i",yaxs="i")#,bg="black")
plot(1,type="n",axes=FALSE,ylim=c(1,100),xlim=c(0,1),xlab="",ylab="")
resolution<-prod(res(meanDensity.R)/1000)
num<-seq(min(cuts)/resolution,max(cuts)/resolution,0.05)
lab<-DoScale(num,1,100)
segments(0,1:length(col),1,1:length(col),col=col,lwd=5,lend=1)
axis(4,lab,num*100)
box()
graphics.off()
}


colFunc <- colorRampPalette(c("gray80","slateblue","yellow","orange","red","red"))
col <- colFunc(100)
{
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_map2.pdf" )),
      width = 10, height = 8)
  par(mfrow = c(1,1), mar = c(4,4,0,0),bg="white")
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       col = "white", border = "black")
  plot( st_geometry(countries), add = T, lwd = 1, border = "black", col = "gray80")
  plot( meanDensity.R, add = T,
        breaks = cuts, col = col,
        axes = F, box = F, bty = "n", legend = F)
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       add=T, border = "black")
  graphics.off()
  
  
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_legend2.pdf" )),
      width = 1.75, height = 5)
  par(mar=c(1,1,1,4),las=1,cex.axis=1.6,xaxs="i",yaxs="i")#,bg="black")
  plot(1,type="n",axes=FALSE,ylim=c(1,100),xlim=c(0,1),xlab="",ylab="")
  resolution<-prod(res(meanDensity.R)/1000)
  num<-seq(min(cuts)/resolution,max(cuts)/resolution,0.05)
  lab<-DoScale(num,1,100)
  segments(0,1:length(col),1,1:length(col),col=col,lwd=5,lend=1)
  axis(4,lab,num*100)
  box()
  graphics.off()
}


colFunc <- colorRampPalette(c("gray80", "yellow", "orange","red","red"))
col <- colFunc(100)
{
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_map3.pdf" )),
      width = 10, height = 8)
  par(mfrow = c(1,1), mar = c(4,4,0,0),bg="white")
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       col = "white", border = "black")
  plot( st_geometry(countries), add = T, lwd = 1, border = "black", col = "gray80")
  plot( meanDensity.R, add = T,
        breaks = cuts, col = col,
        axes = F, box = F, bty = "n", legend = F)
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       add=T, border = "black")
  graphics.off()
  
  
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_legend3.pdf" )),
      width = 1.75, height = 5)
  par(mar=c(1,1,1,4),las=1,cex.axis=1.6,xaxs="i",yaxs="i")#,bg="black")
  plot(1,type="n",axes=FALSE,ylim=c(1,100),xlim=c(0,1),xlab="",ylab="")
  resolution<-prod(res(meanDensity.R)/1000)
  num<-seq(min(cuts)/resolution,max(cuts)/resolution,0.05)
  lab<-DoScale(num,1,100)
  segments(0,1:length(col),1,1:length(col),col=col,lwd=5,lend=1)
  axis(4,lab,num*100)
  box()
  graphics.off()
}


colFunc <- colorRampPalette(c("gray80", "orange", "yellow", "white"))#,"yellow","orange","red","red"
col <- colFunc(100)
{
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_map4.pdf" )),
      width = 10, height = 8)
  par(mfrow = c(1,1), mar = c(4,4,0,0),bg="white")
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       col = "white", border = "black")
  plot( st_geometry(countries), add = T, lwd = 1, border = "black", col = "gray80")
  plot( meanDensity.R, add = T,
        breaks = cuts, col = col,
        axes = F, box = F, bty = "n", legend = F)
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       add=T, border = "black")
  graphics.off()
  
  
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_legend4.pdf" )),
      width = 1.75, height = 5)
  par(mar=c(1,1,1,4),las=1,cex.axis=1.6,xaxs="i",yaxs="i")#,bg="black")
  plot(1,type="n",axes=FALSE,ylim=c(1,100),xlim=c(0,1),xlab="",ylab="")
  resolution<-prod(res(meanDensity.R)/1000)
  num<-seq(min(cuts)/resolution,max(cuts)/resolution,0.05)
  lab<-DoScale(num,1,100)
  segments(0,1:length(col),1,1:length(col),col=col,lwd=5,lend=1)
  axis(4,lab,num*100)
  box()
  graphics.off()
}


colFunc <- colorRampPalette(c("gray80", rev(hcl.colors(4))))
col <- colFunc(100)
{
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_map5.pdf" )),
      width = 10, height = 8)
  par(mfrow = c(1,1), mar = c(4,4,0,0),bg = "white")
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       col = "white", border = "black")
  plot( st_geometry(countries), add = T, lwd = 1, border = "black", col = "gray80")
  plot( meanDensity.R, add = T,
        breaks = cuts, col = col,
        axes = F, box = F, bty = "n", legend = F)
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       add=T, border = "black")
  graphics.off()
  
  
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_legend5.pdf" )),
      width = 1.75, height = 5)
  par(mar=c(1,1,1,4),las=1,cex.axis=1.6,xaxs="i",yaxs="i")#,bg="black")
  plot(1,type="n",axes=FALSE,ylim=c(1,100),xlim=c(0,1),xlab="",ylab="")
  resolution<-prod(res(meanDensity.R)/1000)
  num<-seq(min(cuts)/resolution,max(cuts)/resolution,0.05)
  lab<-DoScale(num,1,100)
  segments(0,1:length(col),1,1:length(col),col=col,lwd=5,lend=1)
  axis(4,lab,num*100)
  box()
  graphics.off()
}


colFunc <- colorRampPalette(c("gray80", hcl.colors(4)))
col <- colFunc(100)
{
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_map6.pdf" )),
      width = 10, height = 8)
  par(mfrow = c(1,1), mar = c(4,4,0,0),bg = "white")
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       col = "white", border = "black")
  plot( st_geometry(countries), add = T, lwd = 1, border = "black", col = "gray80")
  plot( meanDensity.R, add = T,
        breaks = cuts, col = col,
        axes = F, box = F, bty = "n", legend = F)
  plot(st_geometry(st_union(st_intersection(countries,studyArea))),
       add=T, border = "black")
  graphics.off()
  
  
  pdf(file = file.path(thisDir, "figures", paste0(modelName, "_density_legend6.pdf" )),
      width = 1.75, height = 5)
  par(mar=c(1,1,1,4),las=1,cex.axis=1.6,xaxs="i",yaxs="i")#,bg="black")
  plot(1,type="n",axes=FALSE,ylim=c(1,100),xlim=c(0,1),xlab="",ylab="")
  resolution<-prod(res(meanDensity.R)/1000)
  num<-seq(min(cuts)/resolution,max(cuts)/resolution,0.05)
  lab<-DoScale(num,1,100)
  segments(0,1:length(col),1,1:length(col),col=col,lwd=5,lend=1)
  axis(4,lab,num*100)
  box()
  graphics.off()
}




plot( st_geometry(st_difference(countries,studyArea)),
     add= T, lwd = 2, col = "white",border = "white")
plot( st_geometry(bground),
      add = T, lwd = 2, col = "white",border = "white")
# plot( st_geometry(countries), col="black", border = "white", add = T, lwd = 1)
plot( meanDensity.R,
      breaks = cuts, col = col,
      axes = F, box = F, bty = "n", legend = F)
plot( st_geometry(studyArea), add = T, lwd = 3, border = "black")

mtext( text = paste( "N = ", round(WA_regions$summary["Total",1],1),
                     " [", round(WA_regions$summary["Total",4],1), " ; ",
                     round(WA_regions$summary["Total",5],1), "]", sep = ""),
       side = 1, font = 2, cex = 1.5)


##---- Plot Alpine region density raster
par(mfrow = c(1,1), mar = c(4,4,0,0), bg = "black")

alpDensity.R <- habitat.r
alpDensity.R[] <- WA_alps$MeanCell
alpDensity.R[is.na(alps.r[])] <- NA

plot( habitat$polygon, border = "black")
plot( st_geometry(countries), col="black", border = "black", add = T, lwd = 1)
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






## ------     2.1.2. DENSITY EFFECT PLOTS ------
par(mfrow = c(1,1), mar = c(6,6,0,0), bg = "black")
covNames <- names(nimData$hab.covs)
cols1 <- met.brewer(name="Isfahan1",n=length(covNames),type="discrete")
cols2 <- hcl.colors(length(covNames))
cols <- c(cols1[1:3], cols2[2])
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
  plot( x = pred.hab.covs[ ,b], col.lab = "white",
        y = quant.int[2, ], cex.lab = 2,
        type = "n", ylim = c(0, maxD), xlim = range(pred.hab.covs[ ,b]),
        ylab = "Density", xlab = covNames[b], axes = FALSE)
  minCov <- min(st_drop_geometry(habitat$grid[ ,covNames[b]]))
  maxCov <- max(st_drop_geometry(habitat$grid[ ,covNames[b]]))
  xLabels <- round(seq(minCov, maxCov, length.out = 10),2)
  axis(1, col = "white", col.axis = "white",
       at = round(seq(min(pred.hab.covs[ ,b]), max(pred.hab.covs[ ,b]), length.out = 10),3),
       labels = xLabels, cex = 2,
       tck = 0.01, las = 1, hadj = 0.5)
  axis(2, at = seq(0,maxD,length.out = 6), col = "white", col.axis = "white",
       labels = seq(0,maxD,length.out = 6),
       tck = 0.01, las = 1, cex = 2)
  
  polygon(x = c(pred.hab.covs[ ,b],rev(pred.hab.covs[ ,b])),
          y = c(quant.int[1, ],rev(quant.int[3, ])), border = F,
          col = adjustcolor(cols[b], alpha.f = 0.5))
  
  points( x = pred.hab.covs[ ,b],
          y = quant.int[2,],
          lwd = 4, type = "l", col = cols[b])
  
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
par(col.axis = "white", col.lab="white", col.text = "white", bg = "black")

plot(detectors$grid[ ,c("p0_female_alpha", 
                        "p0_male_alpha",
                        "p0_female_pup",
                        "p0_male_pup",
                        "p0_female_other",
                        "p0_male_other")],
     key.pos = 4,
     breaks = c(0,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1)) 



## ------     2.2.2. DETECTION EFFECT PLOTS ------
par(mfrow = c(3,2), mar = c(8,8,0,4))
par(mfrow = c(1,1), mar = c(5,5,0,0), bg = "black",col.lab = "white", col.axis="white")

covNames <- names(nimData$det.covs)
cols1 <- met.brewer(name="Isfahan1",n=length(covNames),type="discrete")
cols2 <- hcl.colors(length(covNames))
cols <- c(cols2[5], cols1[2:4], cols2[2])
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
  axis(1,col = "white",
       at = round(seq(min(pred.det.covs[ ,b]), max(pred.det.covs[ ,b]), length.out = 10),3),
       labels = xLabels, cex = 3,
       tck = 0.01, las = 1, hadj = 0.5)
  axis(2, at = seq(0,maxp0,length.out = 6),col = "white",
       labels = seq(0,maxp0,length.out = 6),
       tck = 0.01, las = 1, cex = 3)
  
  polygon(x = c(pred.det.covs[ ,b],rev(pred.det.covs[ ,b])),
          y = c(quant.int[1, ],rev(quant.int[3, ])), border = F,
          col = adjustcolor(cols[b], alpha.f = 0.5))
  
  points( x = pred.det.covs[ ,b],
          y = quant.int[2,],
          lwd = 5, type = "l", col = cols[b])
  
}


## ------   2.3. SIGMA / p0 ------
par(mfrow = c(1,2))#, bg = "black", col.axis="white", col.text = "white")
#mar = c(15,25,15,25)
myCols <- matrix(met.brewer(name="Isfahan1",n=8,type="discrete")[c(2:4,7:5)],
                 nrow = 3, ncol=2, byrow = F)
valNoise <- c(-0.2,0,0.2)
plot(1,1,
     xlim = c(0.5,2.5),
     ylim = c(0,10),
     type = "n", xaxt = "n", main = "",
     ylab = expression(paste(sigma, " (km)")),
     col.lab = "black",
     xlab = "", axes = F)
axis(1, at = 1:2, labels = sex, col = "black")
axis(2, at = seq(0,10,length.out = 6), col ="black",
     labels = seq(0,10,length.out = 6)*5)
for(s in 1:length(sex)){
  for(ss in 1:length(status)){
    temp <- res$sims.list$sigma[ ,ss,s]
    
    plot.violins2( dat.list = list(temp),
                   x = s + valNoise[ss],
                   at = s + valNoise[ss],
                   add = T,
                   col = myCols[ss,s],
                   alpha = 1,
                   violin.width = 0.2,
                   border.col = myCols[ss,s],
                   scale.width = F)
  }#status
}#sex
# legend(title = "Females",
#        legend = c("Alpha","Pup","Other"),
#        x = 0.5, y = 6.5, fill = myCols[1:3])
# legend(title = "Males",
#        legend = c("Alpha","Pup","Other"),
#        x = 0.9, y = 6.5, fill = myCols[4:6])


plot(1,1,
     xlim = c(0.5,2.5),
     ylim = c(0,0.03),
     type = "n", xaxt = "n", main = "",
     ylab = expression(p[0]),
     col.lab = "black",
     xlab = "", axes = F)
axis(1, at = 1:2, labels = sex, col ="black")
axis(2, at = seq(0,0.03,length.out = 7), col ="black", labels = seq(0,0.03,length.out = 7))
for(s in 1:length(sex)){
  for(ss in 1:length(status)){
    temp <- res$sims.list$p0[ ,ss,s]
    
    plot.violins2( dat.list = list(temp),
                   x = s + valNoise[ss],
                   at = s + valNoise[ss],
                   add = T,
                   col = myCols[ss,s],
                   alpha = 0.9,
                   violin.width = 0.2,
                   border.col = myCols[ss,s],
                   scale.width = F)
  }#status
}#sex
legend(title = "Females",text.col = "black",
       legend = c("Alpha","Pup","Other"),
       x = 0.5, y = 0.03, fill = myCols[1:3])
legend(title = "Males", text.col = "black",
       legend = c("Alpha","Pup","Other"),
       x = 0.9, y = 0.03, fill = myCols[4:6])






## ------   2.3. theta / rho ------
##---- Create empty table
propTable  <-  matrix(NA, nrow = 4, ncol = 3)
colnames(propTable) <- c("female","male","Total")
rownames(propTable) <- c("alpha","pup","other", "Total")
for(ss in 1:length(status)){
    propFemale <- res$sims.list$theta[ ,ss,1] * (1- res$sims.list$rho)
    propMale <- res$sims.list$theta[ ,ss,2] * (res$sims.list$rho)
    propTable[ss,1] <- paste0(
      round( mean(propFemale),2),
      " (", round(quantile(propFemale,0.025), 2),
      "-", round(quantile(propFemale,0.975), 2), ")")
    propTable[ss,2] <- paste0(
      round( mean(propMale),2),
      " (", round(quantile(propMale,0.025), 2),
      "-", round(quantile(propMale,0.975), 2), ")")
    propTable[ss,3] <- paste0(
      round( mean(propFemale+propMale),2),
      " (", round(quantile(propFemale+propMale,0.025), 2),
      "-", round(quantile(propFemale+propMale,0.975), 2), ")")
  }#status

propTable[4,1] <- paste0(
  round( mean(1-res$sims.list$rho),2),
  " (", round(quantile(1-res$sims.list$rho,0.025), 2),
  "-", round(quantile(1-res$sims.list$rho,0.975), 2), ")")

propTable[4,2]<- paste0(
  round( res$mean$rho,2),
  " (", round(res$q2.5$rho, 2),
  "-", round(res$q97.5$rho, 2), ")")





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
# write.csv( popCompo,
#           file = file.path(thisDir, paste0(modelName, "_PopCompo.csv")))



## ------     2.4.2. ABUNDANCES PER REGION/SEX ------
sex <- c("F","M")
status <- c("Alpha","Pups","Others")

##---- Create empty table
abundanceTable  <-  matrix(NA, nrow = 10, ncol = 12)
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
grid.table(abundanceTable[ ,10:12])
mtext( text = "Alpine wolf abundance estimates", 
        side = 3, outer = T, line = -10, font = 2)



##----- Export as .csv
write.csv( abundanceTable[ ,10:12],
          file = file.path(thisDir, paste0(modelName, "_TableAbundance3.csv")))




## ------   2.5. NGS DATA ------
##---- Genotyped NGS samples
ngs <- read.csv(file.path(dataDir,"DNA/ngs21032022.csv"))
ngs$Date <- parse_date_time(ngs$Date, orders = c('dmy','mdy','ymd'))
ngs$Year <- as.numeric(format(ngs$Date,"%Y"))
ngs$Month <- as.numeric(format(ngs$Date,"%m"))
ngs <- ngs[ngs$Year > 2019, ]
ngs <- ngs[ngs$Dead.recovery == "", ]
dim(ngs)

##---- Number of detections per individual
numDetsPerId <- table(ngs$Genotype.ID)
##---- Number of individuals detected
length(numDetsPerId)
##---- Mean number of detections per ID
mean(numDetsPerId)
sd(numDetsPerId)

##---- Turn ngs into spatial data frame
coordinates(ngs) <- cbind.data.frame(ngs$CoordX, ngs$CoordY)
proj4string(ngs) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
ngs <- st_as_sf(ngs)

##---- Plot ngs 
png(filename = file.path(reportDir, "Figures", paste0(modelName, "_NGS_map.png")),
    width = 30, height = 20, units = "cm", res = 600)
par(mar = c(0,0,0,0))
cols <- met.brewer(name="Isfahan1",n=10,type="continuous")
plot( studyArea, border = "black")
plot( st_geometry(countries), col="gray60", border = "white", add = T, lwd = 1)
plot(studyArea,border = "white", col = "gray40", add=T)
plot(transects, add=T, col = "red", lwd = 1)
plot(st_geometry(ngs[ngs$Sex == "", ]), add = T, pch = 3, col = cols[2])
plot(st_geometry(ngs[ngs$Sex == "M", ]), add = T, pch = 3, col = cols[5])
plot(st_geometry(ngs[ngs$Sex == "F", ]), add = T, pch = 3, col = cols[7])
# mtext( text = paste0( "transects : ", round(sum(st_length(transects))/1000),
#                       "km/n genotyped samples : ", dim(ngs)[1]),
#        side = 1, font = 2, cex = 1.5)
legend("bottomright",legend = c("female","male","NA"),pch = 3,col=cols[c(7,5,2)])
dev.off()

#plot( st_geometry(countries), add = T, lwd = 1, border = "white")
plot( st_geometry(studyArea), add = T, lwd = 3, border = "white")




## -----------------------------------------------------------------------------