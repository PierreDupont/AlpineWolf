
################################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone]
##### ---------------- z[psi,rho,theta]
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------ CLEAN THE WORK ENVIRONMENT ------

.libPaths(new=c(.libPaths(),"C:\\DIV\\R-4.0.5\\library"))


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


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.2.2"
thisDir <- file.path(analysisDir, modelName)

## HABITAT SPECIFICATIONS
habitat = list( resolution = 5000, 
                buffer = 30000)

## DETECTORS SPECIFICATIONS
detectors = list( detResolution = 5000,
                  samplingMonths = c(10:12,1:4))

## NGS DATA SPECIFICATIONS
data = list( sex = c("F","M"),
             status = c("alpha","pup","other"),
             aug.factor = 5) 

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
dir.create(thisDir,recursive = TRUE)
dir.create(file.path(thisDir, "input"),recursive = TRUE)
dir.create(file.path(thisDir, "output"),recursive = TRUE)



## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. STUDY AREA ------
##---- Polygon of Italy and neighbouring countries
countries <- read_sf(file.path(dataDir,"GISData/Italy_borders/Italy_andBorderCountries_splitFrance.shp"))

##--- Study area grid
studyAreaGrid <- read_sf(file.path(dataDir,"GISData/shape_studyarea_ALPS/Studyarea_ALPS_2020_2021.shp"))
studyAreaGrid <- st_transform(x = studyAreaGrid, crs = st_crs(countries))
studyArea <- studyAreaGrid %>%
  st_snap(x = ., y = ., tolerance = 0.0001) %>%
  st_union() 

##---- Plot check
plot(st_geometry(studyArea))
plot(st_geometry(countries), col = "gray80",add=T)
plot(studyAreaGrid["SCR"], add = T)









##---- All collected NGS samples

allSamples <- read_sf(file.path(dataDir,"GISData/scats/merged_scats.shp"))

plot(st_geometry(studyArea))

plot(st_geometry(countries), col = "gray80",add=T)

plot(allSamples, add = T,pch=3)



##---- Genotyped NGS samples

ngs <- read.csv(file.path(dataDir,"DNA/ngs.csv"))


ngs.sf <- st_as_sf(x = ngs,                         
               coords = c("CoordX","CoordY"),
               crs = st_crs(allSamples))

plot(st_geometry(ngs.sf[ngs.sf$Type.of.sample%in%c("Scat"),]),add=TRUE,pch=19,col="red",cex=0.6)


##---- POINT PROCESS DENSITY PLOTS





plot(r2)

r<-raster(allSamples,resolution = 50000) 
plot(r)

library(spatstat)

habOwin <- as.owin(as.vector(extent(r)))

pts <- st_coordinates(allSamples)

p <- ppp(pts[,1], pts[,2], window=habOwin)

adjust<-0.5
allsamples.ds <- raster(density(p,adjust=adjust))


plot(allsamples.ds)
#plot(density(log(allsamples.ds[]),na.rm=TRUE))


pts <- st_coordinates(ngs.sf)

p <- ppp(pts[,1], pts[,2], window=habOwin)


genotyped.ds <- raster(density(p,adjust=adjust))

plot(genotyped.ds)
#plot(density(log(genotyped.ds[]),na.rm=TRUE))


ratio<-genotyped.ds/allsamples.ds

ratio[allsamples.ds<quantile(allsamples.ds[],0.3)]<-NA

plot(ratio,axes=FALSE,box=FALSE)

plot(st_geometry(countries), border="black",add=T)
plot(st_geometry(studyArea),add=TRUE,lwd=2)
plot(allSamples, add = T,pch=19,col="grey",cex=0.6)
plot(st_geometry(ngs.sf[ngs.sf$Type.of.sample%in%c("Scat"),]),add=TRUE,pch=19,col="black",cex=0.6)


