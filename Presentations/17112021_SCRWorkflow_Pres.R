####################################################
##### ------- Bavaria National Forest: ------- #####
##### --- Output Processing & Presentation --- #####
####################################################
rm(list=ls())

## ------ LOAD LIBRARIES -----
library(rgdal)
library(raster)
library(rgeos)
library(RANN)
library(R.utils)
library(RColorBrewer)
library(sf)
library(fasterize)
library(nimble)
library(nimbleSCR)
library(coda)
library(sp)
library(plot3D)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(xtable)
library(abind)
library(coda)
library(basicMCMCplots)
library(gridExtra)


## ------ SET YOUR WORKING DIRECTORIES -----
dir.dropboxALL <- "C:/Users/pidu/Dropbox (AQEG)/AQEG Team Folder/"


## ------ SOURCE THE REQUIRED FUNCTIONS ------
source("C:/My_documents/RovQuant/Temp/PD/myWorkingDirectories.R")             
sourceDirectory(dir.function, modifiedOnly = FALSE)
sourceDirectory(dir.function.nimble, modifiedOnly = FALSE)

## WORKING DIRECTORY & MODEL NAME
WD = "C:/Users/pidu/Dropbox (AQEG)/AQEG Team Folder/BavRedDeer/ANALYSES/RED DEER"

## OUTPUT DIRECTORY 
outFolder = "Results_AC[exp(REG+DEM.R)]_Z[.]_Y[admin+tracks]"

## -----------------------------------------------------------------------------
## ------ I. LOAD AND PREPARE DATA ------
load(file.path( myVars$WD,
                myVars$outFolder,
                "processedObjects.RData"))
poly <- readOGR(file.path(dir.dropboxALL, "BavRedDeer/cleaned/NR_Cleaned_WGS84Z33N.shp"))
grid <- readOGR(file.path(dir.dropboxALL, "BavRedDeer/Data/Search_grids_NR.shp"))
grid <- spTransform(grid, crs(poly))

load(file = file.path(dir.dropboxALL, "BavRedDeer/EffortData.RData")) ## "Tracks
tracks <- gSimplify(Tracks,200)
tracks <- tracks[!is.na(over(tracks,poly)[,1]), ]
plot(poly, col = "gray60", border = F)
plot(grid, add = T)
plot(tracks, add = T, col = "red")


## DETECTORS
detectors1000 <- raster(extent(tracks), resolution = 1000)
detectors1000 <- rasterize(tracks, detectors1000, fun = "length")

detectors750 <- raster(extent(tracks), resolution = 750)
detectors750 <- rasterize(tracks, detectors750, fun = "length")

detectors500 <- raster(extent(tracks), resolution = 500)
detectors500 <- rasterize(tracks, detectors500, fun = "length")


## HABITAT 
habitat <- gBuffer(spgeom = aggregate(rasterToPolygons(detectors500)),width = 2000)
hab.r <- raster(extent(habitat), resolution = 1000)
hab.r <- rasterize(habitat, hab.r)
elev.r <- hab.r
elev.r[ ] <- raster::extract(x = DEM.r,
                             y = SpatialPoints(coords = coordinates(elev.r),
                                               proj4string = crs(elev.r)))


elev.r[is.na(hab.r)] <- NA
## -----------------------------------------------------------------------------





## DETECTORS
png(file= "C:/Users/pidu/StudyArea.png", width = 10,height = 10,res = 600,units = "in")
par(bg=NA)
plot(poly, col = "gray60", border = F)
dev.off()

png(file= "C:/Users/pidu/StudyGrid.png", width = 10,height = 10,res = 600,units = "in")
par(bg=NA)
plot(poly, col = "gray60", border = F)
plot(grid, border = rgb(red = 0,green = 154/283,blue = 129/283), add=T)
dev.off()

png(file= "C:/Users/pidu/StudyTracks.png", width = 10,height = 10,res = 600,units = "in")
par(bg = NA)
plot(poly, col = "gray60", border = F)
plot(grid, border = rgb(red = 0,green = 154/283,blue = 129/283), add=T)
plot(tracks, add = T,col ="red")
dev.off()

png(file= "C:/Users/pidu/Detectors1000.png", width = 10,height = 10,res = 600,units = "in")
par(bg=NA)
plot(poly, col = "gray60", border = F)
plot(rasterToPolygons(detectors1000),
     border = rgb(red = 0,green = 154/283,blue = 129/283),
     col = adjustcolor("gray40",0.4), add = T)
dev.off()

png(file= "C:/Users/pidu/Detectors750.png", width = 10,height = 10,res = 600,units = "in")
par(bg=NA)
plot(poly, col = "gray60", border = F)
plot(rasterToPolygons(detectors750),
     border = rgb(red = 0,green = 154/283,blue = 129/283),
     col = adjustcolor("gray40",0.4), add = T)
dev.off()

png(file= "C:/Users/pidu/Detectors500.png", width = 10,height = 10,res = 600,units = "in")
par(bg=NA)
plot(poly, col = "gray60", border = F)
plot(rasterToPolygons(detectors500),
     border = rgb(red = 0,green = 154/283,blue = 129/283),
     col = adjustcolor("gray40",0.4), add = T)
dev.off()


png(file= "C:/Users/pidu/Detectors.png", width = 10,height = 10,res = 600,units = "in")
par(bg=NA)
plot(poly, col = "gray60", border = F)
image(detectors500, border = rgb(red = 0,green = 154/283,blue = 129/283), add=T)
plot(rasterToPolygons(detectors500),
     border = rgb(red = 0,green = 154/283,blue = 129/283),
     add = T)
dev.off()



## HABITAT
png(file= "C:/Users/pidu/Habitat1.png", width = 10,height = 10,res = 600, units = "in")
par(bg=NA)
plot(habitat, col = adjustcolor("forestgreen",0), border = F)
plot(poly, border = rgb(red = 0,green = 154/283,blue = 129/283),add=T)
plot(rasterToPolygons(detectors500), add=T, border = "white")
dev.off()

png(file= "C:/Users/pidu/Habitat2.png", width = 10,height = 10,res = 600, units = "in")
par(bg=NA)
plot(habitat, col = adjustcolor("forestgreen",0.6),
     border = rgb(red = 0,green = 154/283,blue = 129/283))
plot(poly, border = rgb(red = 0,green = 154/283,blue = 129/283),add=T)
plot(rasterToPolygons(detectors500),add=T,border = "white")
dev.off()

png(file= "C:/Users/pidu/Habitat1000.png", width = 10,height = 10,res = 600,units = "in")
par(bg=NA)
plot(rasterToPolygons(hab.r), col = adjustcolor("forestgreen",0.6),
     border = rgb(red = 0,green = 154/283,blue = 129/283))
plot(poly, border = rgb(red = 0,green = 154/283,blue = 129/283),add=T)
plot(rasterToPolygons(detectors500),add=T,border = "white")
dev.off()

png(file= "C:/Users/pidu/Elevation1000.png", width = 10,height = 10,res = 600,units = "in")
par(bg = NA)
plot(rasterToPolygons(hab.r))
plot(elev.r, border = rgb(red = 0,green = 154/283,blue = 129/283),add=T)
plot(rasterToPolygons(elev.r), add=T,
     border = rgb(red = 0,green = 154/283,blue = 129/283))
plot(rasterToPolygons(detectors500),add=T,border = "white")
plot(poly, add=T)
dev.off()

