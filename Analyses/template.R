######################################################
##### ------------- ALPINE WOLD SCR ------------ #####
##### ----- D[...] ------- #####
##### ----- p0[...] ------ #####
##### ----- sigma[] ------ #####
###################################################
## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())

## ------ IMPORT REQUIRED LIBRARIES ------
library(rgdal)
library(raster)
library(sp)
library(coda)
library(nimble)
library(nimbleSCR)
library(spdep)
library(rgeos)
library(maptools)
library(stringr)
library(abind)
library(R.utils)
library(adehabitatHR)
library(sf)
library(fasterize)

## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")

## ------ SOURCE CUSTOM FUNCTIONS (IF NECESSARY ------



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## Model name
modelName = "AW0.0"
WD = analysisDir

## HABITAT SPECIFICATIONS
habitat = list( polygon =  c("SWE","NOR"),
                resolution = 20000, 
                buffer = 60000)

## NGS DATA SPECIFICATIONS
dna = list( years = 2012:2020,   
            species = c("Jerv"),               ## "Ulv","Jerv","Bjorn"
            sex = c("Hunn"),                   ## "Hann","Hunn","Ukjent" 
            samplingMonths = list(12,1:6)) ## list(10:12,1:4), list(1:XXX), list(XX:XX,YY:YY)

## DETECTORS SPECIFICATIONS
detectors = list( detSubResolution = 2000,
                  detResolution = 10000,
                  detDeadResolution = 15000)

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(is.null(WD))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
if(!dir.exists(file.path(WD, modelName))){dir.create(file.path(WD, modelName))}

## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
### ==== 1. HABITAT DATA ====
### ====    1.1.LOAD RAW SHAPEFILES ====
## POLYGONS OF THE REGION
GLOBALMAP <- readOGR(paste(dir.dropbox,"/DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp",sep="")) ## Map of Scandinavia (including Finland & parts of Russia)
GLOBALMAP <- GLOBALMAP[GLOBALMAP$area > 80000000, ]
GLOBALMAP <- crop(GLOBALMAP, extent(c(-70000,1200000,5100000,8080000)))

## POLYGONS OF SWEDEN & NORWAY
COUNTRIES <- GLOBALMAP[GLOBALMAP$ISO %in% c("SWE","NOR"), ]
COUNTRIES <- aggregate(x = COUNTRIES, by = "ISO")

## POLYGONS OF COMMUNES IN SWEDEN & NORWAY
COMMUNES_NOR <- readOGR(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp", sep = ""))   ## Communal map of Norway
COMMUNES_SWE <- readOGR(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp", sep = ""))    ## Communal map of Sweden
COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)
## POLYGONS OF COUNTIES IN SWEDEN & NORWAY
COUNTIES <- aggregate(x = COMMUNES, by = "NAME_1")

## AGGREGATE COUNTIES (OPTIONAL)
COUNTIES_AGGREGATED <- COUNTIES
COUNTIES_AGGREGATED <- gSimplify(COUNTIES_AGGREGATED,tol=500, topologyPreserve = TRUE)
COUNTIES_AGGREGATED$id <- 1:length(COUNTIES_AGGREGATED)
#[CM] adjust Counties aggregation 
COUNTIES_AGGREGATED$id[c(2:3,6,24,40,13)] <- 3
COUNTIES_AGGREGATED$id[c(35)] <- 14
COUNTIES_AGGREGATED$id[c(14,36)] <- 15
COUNTIES_AGGREGATED$id[c(15,16,1,23)] <- 15
COUNTIES_AGGREGATED$id[c(27,37,33,30)] <- 15
COUNTIES_AGGREGATED$id[c(38,34,7,9)] <- 15
COUNTIES_AGGREGATED$id[c(26,19)] <- 2
COUNTIES_AGGREGATED$id[c(22,18,12)] <- 3
COUNTIES_AGGREGATED$id[c(29,2,31,25,4,39)] <- 3#[CM]

COUNTIES_AGGREGATED <- aggregate(x = COUNTIES_AGGREGATED, by = "id")

plot(gSimplify(COUNTIES_AGGREGATED,tol=500), border=grey(0.5), col="grey")
text(COUNTIES_AGGREGATED, labels = COUNTIES_AGGREGATED$id, col = "black")  

### ====    1.3.SAVE SHAPEFILES OBJECTS FOR FASTER RUNS ====
# save( GLOBALMAP,
#       COUNTRIES,
#       COUNTIES,
#       COUNTIES_AGGREGATED,
#       COMMUNES,
#       file = file.path(myVars$WD, myVars$modelName, "HABITAT.RData"))
# load(file.path(myVars$WD, myVars$modelName, "HABITAT.RData"))

### ====    1.4.CREATE STUDY AREA POLYGON ====
## CREATE STUDY AREA POLYGON BASED ON COUNTY & COMMUNES IDs
if(!is.null(myVars$HABITAT$countyNames)){
   myStudyArea <- COMMUNES[COMMUNES$NAME_1 %in% myVars$HABITAT$countyNames, ]
   myStudyArea1 <- COMMUNES[COMMUNES$ISO == "NOR" & COMMUNES$ID_2 %in% myVars$HABITAT$communeNames, ]
   myStudyArea <- aggregate(rbind(myStudyArea, myStudyArea1)) 
   ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
   myBufferedArea <- gBuffer(spgeom = myStudyArea, width = myVars$HABITAT$habBuffer)
   myBufferedArea <- gIntersection(myBufferedArea, GLOBALMAP)
}

## CREATE STUDY AREA POLYGON BASED ON x AND y EXTENTS
if(!is.null(myVars$HABITAT$x.extent)){
   myStudyArea <- crop(COUNTRIES, extent(myVars$HABITAT$x.extent, myVars$HABITAT$y.extent))
   ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
   myBufferedArea <- gBuffer(spgeom = myStudyArea, width = myVars$HABITAT$habBuffer)
   myBufferedArea <- gIntersection(myBufferedArea, GLOBALMAP)
   myBufferedArea <- crop(myBufferedArea, extent(myVars$HABITAT$x.extent, myVars$HABITAT$y.extent))
}

## CREATE STUDY AREA POLYGON BASED ON COUNTRY NAMES
if(!is.null(myVars$HABITAT$countries)){
   myStudyArea <- COUNTRIES[COUNTRIES$ISO %in% myVars$HABITAT$countries, ]
   ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
   myBufferedArea <- gBuffer(spgeom = myStudyArea, width = myVars$HABITAT$habBuffer)
   myBufferedArea <- gIntersection(myBufferedArea, GLOBALMAP)
}

## CREATE STUDY AREA POLYGON BASED ON POLYGON COORDINATES
if(!is.null(myVars$HABITAT$coord.x)){
   myStudyPoly <- MakePolygon( coord.x = myVars$HABITAT$coord.x,
                               coord.y = myVars$HABITAT$coord.y)
   myStudyArea <- gIntersection(myStudyPoly, COUNTRIES)
   ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
   myBufferedArea <- gBuffer(spgeom = myStudyArea, width = myVars$HABITAT$habBuffer)
   myBufferedArea <- gIntersection(myBufferedArea, GLOBALMAP)
   myBufferedArea <- gIntersection(myBufferedArea, myStudyPoly)
}

## PLOT CHECK
if(myVars$plot.check){
   par(mfrow = c(1,1))
   plot(COUNTRIES)
   plot(myBufferedArea, add = TRUE, col = rgb(0.72,0.14,0.14,0.3))
   plot(myStudyArea, add = TRUE, col ="red")
}

myStudyArea.poly <- myStudyArea #[CM] to get the tracks to work
### ==== 2. NGS DATA ====
### ====    2.1.LOAD ROVBASE FILES ====
# DNA <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20191008/dna_wolf.csv",sep=""))## NGS data from RovBase
# DEAD <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20191008/dead_carnivores.csv",sep="")) ## Dead Recoveries from RovBase
DNA <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20211029/dna_wolverine.csv", sep=""))## NGS data from RovBase#[CM update to 20190627]
DEAD <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20211029/dead_carnivores.csv",sep="")) ## Dead Recoveries from RovBase#[CM update to 20190627]
SUSPECT_NGS_SAMPLES <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20211029/Remove ngs samples list wolverine.csv", sep = "")) ## DNA samples to be removed from Henrik
SUSPECT_DeadRecoSAMPLES <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20211029/Remove dead recoveries list wolverine.csv", sep = "")) ## DNA samples to be removed from Henrik

DEN <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20211029/DEN_COUNTS_2009_2021_fromHB.csv", sep = ""))

### ====    2.2.TRANSLATE SCANDINAVIAN CHARACTERS ====
colnames(DNA) <- translateForeignCharacters(dat=colnames(DNA), dir.translation = dir.analysis )
colnames(DEAD) <- translateForeignCharacters(dat=colnames(DEAD), dir.translation = dir.analysis )
colnames(DEN) <- translateForeignCharacters(dat = colnames(DEN), dir.translation = dir.analysis)

### ====    2.2.SAVE NGS OBJECTS FOR FASTER RUNS ====
# save(DNA, DEAD, INDIVIDUAL_ID, file = file.path(myVars$WD,"NGS_DATA.RData"))
# load(file.path(myVars$WD,"NGS_DATA.RData"))

## ==== 3. SEARCH EFFORT DATA ====
### ====    3.1.GPS SEARCH TRACKS ====
# LOAD GPS SEARCH TRACKS FROM ROVBASE
# TRACKS_SINGLE <- readOGR(paste(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20211029/AktivitetsloggRovQuant20211029/dbo.XX_eksport_rovquant_aktivitetslogg_alle_spor_20211029_spatial_linestring_date.shp", sep = ""))
# TRACKS_MULTI <- readOGR(paste(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20211029/AktivitetsloggRovQuant20211029/dbo.XX_eksport_rovquant_aktivitetslogg_alle_spor_20211029_spatial_multilinestring_date.shp", sep = ""))
# ## COMBINE ALL TRACKS
# ALL_TRACKS <- rbind(TRACKS_SINGLE, TRACKS_MULTI)
# ## REMOVE HELICOPTER TRACKS
# ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Helkptr == "N", ]
# ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Jerv == "Y", ] #[CM}
# ## SELECT TRACKS YEAR
# TRACKS_YEAR <- TRACKS_YEAR.sp <-list()
# for(t in 1:nYears){
#    ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
#    TRACKS_1 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][1] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[1]], ]
#    TRACKS_2 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][2] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[2]], ]
#    TRACKS <- rbind(TRACKS_1, TRACKS_2)
#    ## SIMPLIFY TRACKS SHAPES
#    TRACKS <- gSimplify(spgeom = TRACKS, tol = 500)
#    ## SUBSET TRACKS TO THE STUDY AREA
#    TRACKS <- intersect(TRACKS, myStudyArea.poly)
#    ## NAME TRACKS
#    TRACKS$ID <- 1:length(TRACKS)
#    TRACKS_YEAR[[t]] <- TRACKS
#    ## TRANSFORM TRACKS INTO SPATIAL POINTS
#    TRACKS_YEAR.sp[[t]] <- SampleTrack(x = TRACKS_YEAR[[t]], d = 200, poly = NULL)
# }#t
# 
# ### ====    3.2.DISTANCE TO ROADS ====
# ## LOAD MAP OF DISTANCES TO ROADS (1km resolution)
# DistAllRoads <- raster(paste(dir.dropbox,"/DATA/GISData/Roads/MinDistAllRoads1km.tif", sep=""))
# #[CM] ACTIVATE!
# r <- fasterize(st_as_sf(myStudyArea.poly), DistAllRoads)
# r[!is.na(r)] <- DistAllRoads[!is.na(r)]
# DistAllRoads <- r
# DistAllRoads <- crop(DistAllRoads,myStudyArea.poly)#[CM]extent(myStudyArea.poly))
# ## PLOT CHECK
# if(myVars$plot.check){
#    plot(DistAllRoads)
#    plot(myStudyArea.poly,add=T)
# }
# 
# ### ====    3.3.DAYS OF SNOW ====
# ## SEASONAL MAPS (CREATED IN TEMP/CM/GIS/snowMODIS)
#  SNOW <- stack(paste(dir.dropbox,"/DATA/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2008_2021.tif", sep=""))
# ## RENAME THE LAYERS
#  names(SNOW) <- paste(2008:2020,(2008:2020)+1, sep="_")
# ## SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
#  SNOW <- SNOW[[paste("X", years, "_", years+1, sep="")]]
#  SNOW <- raster::crop(SNOW, c(0,40,55,75))
# 
# ### ====    3.4.SAVE SEARCH EFFORT OBJECTS FOR FASTER RUNS ====
# save(TRACKS_YEAR, TRACKS_YEAR.sp, SNOW, DistAllRoads, file = file.path(myVars$WD, "TRACKS.RData"))
load(file.path(myVars$WD, "TRACKS.RData"))



## -----------------------------------------------------------------------------
## ------ II. PREPARE SCR DATA ------
### ==== 1.CLEAN AND FILTER NGS DATA ====
## Remove DEAD entries from the DNA data [HB]
DNA <- DNA[substr(DNA$RovbaseID..Proeve.,1,1) != "M", ]

## Remove un-verified dead recoveries [HB] 
## ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
DEAD <- DEAD[!grepl(pattern = "belastet kvote", x = as.character(DEAD$Doedsarsak)), ]




### ====    1.1.CLEAN NGS & DEAD RECOVERY DATA ====
myCleanedData.sp <- CleanDataNew1( dna_samples = DNA
                                   ,
                                   dead_recoveries = DEAD
                                   ,
                                   species_id = myVars$DATA$species
                                   ,
                                   country_polygon = COUNTRIES
                                   ,
                                   threshold_month = unlist(myVars$DATA$samplingMonths)[1]
                                   ,
                                   keep_dead = T
                                   ,
                                   age.label.lookup = age.lookup.table
)


## PLOT CHECK
if(myVars$plot.check){
   plot(COUNTRIES)
   plot(myStudyArea.poly, add = T, col ="red")
   plot(myCleanedData.sp, add=TRUE,pch=19,cex=0.2, col="white")
}


### ====    1.2. FILTER DATA FOR SEX myFullData.sp====
myFullData.sp <- FilterData( myData = myCleanedData.sp,
                             poly = myStudyArea.poly,
                             dead.recovery = T,
                             sex = myVars$DATA$sex,
                             setSex = T)

if(myVars$plot.check){
   plot(myFullData.sp$alive, add=TRUE,pch=19,cex=0.2, col="lightblue")
}

# source("C:/My_documents/RovQuant/Temp/PD/Source/CleanData.R")
# source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/PD/Source/CleanData.R")
# 
# cleanData <- CleanData( dna_samples = DNA,
#                         dead_recoveries = DEAD,
#                         species_id = myVars$DATA$species,
#                         country_polygon = COUNTRIES,
#                         threshold_month = unlist(myVars$DATA$samplingMonths)[1],
#                         setSex = TRUE,
#                         keep_dead = TRUE,
#                         reset_month = FALSE,             
#                         keep_unidentified_id = FALSE,
#                         age.label.lookup = NULL)

## REMOVE SUSPECT SAMPLES ACCORDING TO HENRIK
myFullData.sp$alive$DNAID <- as.character(myFullData.sp$alive$DNAID)
myFullData.sp$dead.recovery$DNAID <- as.character(myFullData.sp$dead.recovery$DNAID)

myFullData.sp$alive <- myFullData.sp$alive[!(myFullData.sp$alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[!(myFullData.sp$dead.recovery$RovBaseId %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]



##EXPORT THE DATA 
if(myVars$DATA$sex=="Hann"){
   assign("myFullData.spM", myFullData.sp)
}else{
   assign("myFullData.spF", myFullData.sp)
}



## Remove individuals that died twice# [CM] TO BE CHECKED BECAUSE "length(IdDoubleDead) < 0" and so it was desactivated
IdDoubleDead <- myFullData.sp$dead.recovery$Id[duplicated(myFullData.sp$dead.recovery$Id)]
if(length(IdDoubleDead) > 0){
   duplicatedDeath <- NULL
   for(i in IdDoubleDead){
      tmp  <- which(myFullData.sp$dead.recovery$Id == i & is.na(myFullData.sp$dead.recovery$DeathCause_2))
      if(length(tmp)==0){tmp  <- which(myFullData.sp$dead.recovery$Id == i)[-1]}
      duplicatedDeath <- c(duplicatedDeath, tmp)
   }#i  
   myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-duplicatedDeath, ]
}#if
myFullData.sp$dead.recovery@data <- droplevels(myFullData.sp$dead.recovery@data)


unique(myFullData.sp$dead.recovery$DeathCause)
unique(myFullData.sp$dead.recovery$DeathCause_2)


## Remove pups killed before recruitment based on weight (cf. Henrik)
## 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and recovered dead between March and November
sum(myFullData.sp$dead.recovery$Doedt.individ..Unge %in% "Ja" &
       myFullData.sp$dead.recovery$Month > 2 &
       myFullData.sp$dead.recovery$Month < 12)
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$Doedt.individ..Unge %in% "Ja" &
                                                                     myFullData.sp$dead.recovery$Month > 2 &
                                                                     myFullData.sp$dead.recovery$Month < 12),]


# 2) remove individuals that have a weight >0 and <4 between March and November 
# format the weight correctly 
myFullData.sp$dead.recovery$Helvekt <- as.character(myFullData.sp$dead.recovery$Helvekt)
myFullData.sp$dead.recovery$Slaktevekt <- as.character(myFullData.sp$dead.recovery$Slaktevekt)

#convert to decimals
myFullData.sp$dead.recovery$Helvekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Helvekt))
myFullData.sp$dead.recovery$Slaktevekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Slaktevekt))
# get the two weight columns together. 
myFullData.sp$dead.recovery$weight <- ifelse(!is.na(myFullData.sp$dead.recovery$Helvekt), 
                                             myFullData.sp$dead.recovery$Helvekt, 
                                             myFullData.sp$dead.recovery$Slaktevekt)
#assign negative values to nas to avoid issues 
myFullData.sp$dead.recovery$weight[is.na(myFullData.sp$dead.recovery$weight)] <- -999
# check how many dead reco we remove and remove if more than 0
if(sum(myFullData.sp$dead.recovery$weight > 0 &
       myFullData.sp$dead.recovery$weight < 4 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2)>0){
   myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$weight > 0 &
                                                                        myFullData.sp$dead.recovery$weight < 4 &
                                                                        myFullData.sp$dead.recovery$Month < 12 &
                                                                        myFullData.sp$dead.recovery$Month > 2),]
}
# check how many dead reco with a weight of 0 kg and recovered between march and november
if(sum(myFullData.sp$dead.recovery$Age %in% 0 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2)>0){
   myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Age %in% 0 &
                                  myFullData.sp$dead.recovery$Month < 12 &
                                  myFullData.sp$dead.recovery$Month > 2,  ]
}

## IDENTIFY PUPS BASED ON AGE AND WEIGHT (old version )
# myFullData.sp$dead.recovery$weight <- ifelse(!is.na(myFullData.sp$dead.recovery$Helvekt), myFullData.sp$dead.recovery$Helvekt, myFullData.sp$dead.recovery$Slaktevekt)
# myFullData.sp$dead.recovery$Age[is.na(myFullData.sp$dead.recovery$Age)] <- -999
# myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[which(!(myFullData.sp$dead.recovery$Month %in% c(3:9) & myFullData.sp$dead.recovery$Age == 0)), ]
# myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[which(myFullData.sp$dead.recovery$weight >= 6), ]
# tmp <- myFullData.sp$dead.recovery[which(!(myFullData.sp$dead.recovery$Month %in% c(3:9) & myFullData.sp$dead.recovery$Age == 0)), ]
# table(tmp$Country, tmp$Year)
# table(myFullData.sp$dead.recovery$Country, myFullData.sp$dead.recovery$Year)




### ====    1.2.FILTER NGS & DEAD RECOVERY DATA ====
myFilteredData.sp <- myFullData.sp

## Subset to years of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years, ]
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year %in% years, ]

## Subset to months of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Month %in% unlist(myVars$DATA$samplingMonths), ]
# myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Month %in% unlist(myVars$DATA$samplingMonths), ]

## Remove dead recoveries outside the HABITAT #[CM] DO this after creating myhabitat
# deadData <- deadData[which(!is.na(over(deadData, as(myBufferedArea,"SpatialPolygons")))), ]

## PLOT CHECK
if(myVars$plot.check){
   par(mfrow = c(1,3))
   for(t in 1:nYears){
      ## DEAD RECOVERIES TOTAL
      tempTotal <-  myFilteredData.sp$dead.recovery[ myFilteredData.sp$dead.recovery$Year == years[t] & myFilteredData.sp$dead.recovery$Sex %in% myVars$DATA$sex, ]
      NGS_TabTotal <- table(tempTotal$Country)
      ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x) sum(x>0))
      ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
      # tempIn <- deadData[deadData$Year == years[t], ]
      # NGS_TabIn <- table(tempIn$Country)
      # ID_TabIn <- apply(table(tempIn$Id, tempIn$Country), 2, function(x) sum(x>0))
      ## PLOT NGS SAMPLES
      plot(GLOBALMAP, col="gray80")
      plot(myStudyArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
      plot(myBufferedArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
      points(tempTotal, pch = 21, bg = "darkred")
      # points(tempIn, pch = 21, bg = "darkgreen")
      # ## ADD NUMBER OF NGS samples and IDs per COUNTRY
      graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
      graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
      ## ADD OVERALL NUMBERS
      mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
   #    mtext(text = paste(sum(NGS_TabIn), "Dead Recoveries /", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
   #    mtext(text = paste(sum(NGS_TabTotal)-sum(NGS_TabIn), "Recoveries /", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"), side = 3, line = -1)
    }#t
}

### ==== 2. GENERATE HABITAT ====
### ====        1.2.1.1 TRY TO REDUCE THE AREA OF THE STATE-SPACE BASED ON DETECTIONS ====
## DELINEATE IT BASED ON DETECTIONS AND DEAD RECOVERY
par(mar=c(0,0,0,0))
plot(myStudyArea)
points(myFilteredData.sp$alive, pch=21, bg="red", cex=0.5)
points(myFilteredData.sp$dead.recovery, pch=21, bg="blue", cex=0.5)

## DELINEATE A BUFFER AROUND ALL DEADRECO 
BuffDead <- aggregate(gBuffer(myFilteredData.sp$dead.recovery, byid = TRUE, width = myVars$HABITAT$habBuffer))
plot(aggregate(BuffDead),add=T, border="blue")
## DELINEATE A BUFFER AROUND ALL DETECTIONS 
BuffAlive <- gBuffer(myFilteredData.sp$alive, byid = TRUE, width = myVars$HABITAT$habBuffer*1.4)
myBufferedArea <- aggregate(BuffAlive)
plot(myBufferedArea, border="red", add=T)

## CUT TO SWEDISH AND NORWEGIAN BORDERS
myStudyArea <- raster::intersect(myBufferedArea, myStudyArea.poly[myStudyArea.poly$ISO %in% c("SWE","NOR"), ])
plot(myStudyArea, border="grey", add=T)
## remove holes so samples in lakes are not removed
# BuffAliveNoHoles <- RemoveHolesSp(myBufferedArea)
# plot(BuffAliveNoHoles, border="yellow", add=T)

### ====    2.1. GENERATE HABITAT CHARACTERISTICS ====
## CREATE HABITAT LIST CHARACTERISTICS
# source("C:/My_documents/RovQuant/Temp/PD/Source/MakeHabitat.R")
# source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/PD/Source/MakeHabitat.R")

myHabitat <- MakeHabitat( poly = as(myStudyArea,"SpatialPolygons"),
                          resolution = myVars$HABITAT$habResolution,
                          buffer = myVars$HABITAT$habBuffer,
                          polygon.clip = GLOBALMAP,#[CM]
                          plot.check = FALSE,
                          fasterize = TRUE,
                          CoverToKeepHabitat = 51)#[CM]# 21

proj4string(myHabitat$buffered.habitat.poly) <- CRS(proj4string(myHabitat$habitat.sp))

temp <- rasterToPolygons(myHabitat$habitat.r)
myHabitat$buffered.habitat.poly <- aggregate(temp[myHabitat$habitat.r[] == 1, ])

## RETRIEVE HABITAT WINDOWS BOUNDARIES
lowerHabCoords <- coordinates(myHabitat$habitat.r)[myHabitat$habitat.r[]==1,] - 0.5*myVars$HABITAT$habResolution
upperHabCoords <- coordinates(myHabitat$habitat.r)[myHabitat$habitat.r[]==1,] + 0.5*myVars$HABITAT$habResolution
nHabCells <- dim(lowerHabCoords)[1]

## create habitat grid 
habIDCells.mx <- myHabitat$IDCells.mx 
habIDCells.mx[] <- 0
scaledHabGridCenters <- UTMToGrid(data.sp = myHabitat$habitat.sp,
                                  grid.sp = myHabitat$habitat.sp,plot.check = F)$grid.scaled.xy

scaledHabGridCenters <- scaledHabGridCenters[myHabitat$habitat.r[]==1,]
for(i in 1:nrow(scaledHabGridCenters)){
   habIDCells.mx[trunc(scaledHabGridCenters[i,2])+1,
                 trunc(scaledHabGridCenters[i,1])+1] <- i
}
image(habIDCells.mx)


### ====    2.1.1 SUBSET DETECTIONS BASED ON HABITAT EXTENT ==== #[CM]
## Remove samples outside the STUDY AREA #[CM]
myFilteredData.sp$alive <- myFilteredData.sp$alive[!is.na(over(myFilteredData.sp$alive, as(myStudyArea,"SpatialPolygons"))), ]
myFilteredData.sp$alive@data <- droplevels(myFilteredData.sp$alive@data)

## Remove dead recoveries outside the HABITAT #[CM] 
myHabitat$buffered.habitat.polySp <- as(myHabitat$buffered.habitat.poly,"SpatialPolygons")
proj4string(myHabitat$buffered.habitat.polySp) <- CRS(proj4string(myHabitat$habitat.sp))

myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[which(!is.na(over(myFilteredData.sp$dead.recovery,
                                                                                     myHabitat$buffered.habitat.polySp ))), ]
myFilteredData.sp$dead.recovery@data <- droplevels(myFilteredData.sp$dead.recovery@data)

## PLOT CHECK
if(myVars$plot.check){
   # par(mfrow = c(1,2))#[CM]
   plot(myHabitat$habitat.r)
   plot(myStudyArea, add = T, col = rgb(150/250,150/250,150/250, alpha = 0.75))
   plot(GLOBALMAP, add = T)
   plot(myHabitat$buffered.habitat.poly, add=T)
   points(myFilteredData.sp$alive,pch=21, bg="red", cex=0.5)
   points(myFilteredData.sp$dead.recovery,pch=21, bg="blue", cex=0.5)
}

### ====    2.2. GENERATE HABITAT-LEVEL COVARIATES ====
### ====       2.2.1. DEN COUNTS ====
DEN.sp <- SpatialPoints(DEN[ ,c("UTM33_X" ,"UTM33_Y")], proj4string = CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# DEN.r <- myHabitat$habitat.r
# DEN.r <- rasterize(DEN.sp, DEN.r, fun = 'count')
# DEN.r[is.na(DEN.r)] <- 0
# denCounts <- DEN.r[myHabitat$habitat.r[] == 1]
# if(myVars$plot.check){
#    plot(DEN.r,main="Sum Den-2013:2018")
#    plot(myStudyArea, add = TRUE, border = "black")
# }

# ##[CM]
# pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"Sum Den20132018Kernel.pdf",sep="")))
# H=seq(5000,60000,by=5000)
# for(h in 1:length(H)){
DEN.sp$id  <- rep(1,length(DEN.sp))
DEN.r <- raster(
   estUDm2spixdf(
      kernelUD(DEN.sp, h = 30000,#H[h],#[CM]
               grid = as(myHabitat$habitat.r, 'SpatialPixels'))
   )
)

if(myVars$plot.check){
   plot(DEN.r)#,main=paste("Sum Den-2013:2018Kernel H=", H[h], sep=""))
   plot(myStudyArea, add = TRUE, border = "black")
}
# }
# dev.off()
# yearly Alternative #[CM}
# DEN.sp$id  <- DEN$ar
# denID <- unique(DEN.sp$id)
# DEN.rYear <- list()
# pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"Year Den20132018Kernel.pdf",sep="")))
# for(t in 1:length(denID)){
#    DEN.rYear[[t]] <- raster(
#    estUDm2spixdf(
#       kernelUD(DEN.sp[DEN.sp$id==denID[t],], h = 30000,#[CM]
#                grid = as(myHabitat$habitat.r, 'SpatialPixels'))
#    )
# )
# if(myVars$plot.check){
#    plot(DEN.rYear[[t]],main=paste(years[t]," H=", 30000, sep="") )
#    plot(myStudyArea, add = TRUE, border = "black")
# }
# }
# dev.off()
##EXTRACT COVARIATEs
denCounts <- DEN.r[myHabitat$habitat.r[ ]==1]
denCounts <- round(scale(denCounts),digits = 2)#[CM]

### ==== 3. GENERATE DETECTORS ====
### ====    3.1. GENERATE DETECTORS CHARACTERISTICS ====
## GENERATE NGS DETECTORS BASED ON THE STUDY AREA
myDetectors <- MakeSearchGrid( data = myStudyArea,
                               resolution = myVars$DETECTORS$detResolution,
                               div = (myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution)^2,
                               plot = FALSE,
                               fasterize = TRUE)

## GENERATE DEAD RECOVERY DETECTORS BASED ON THE AVAILABLE HABITAT
myDetectors.dead <- MakeSearchGrid( data = myHabitat$buffered.habitat.poly,
                                    resolution = myVars$DETECTORS$detDeadResolution,
                                    plot = FALSE,
                                    fasterize = TRUE)

## EXTRACT NUMBERS OF DETECTORS
n.detectors <- dim(myDetectors$main.detector.sp)[1]
n.detectors.dead <- dim(myDetectors.dead$main.detector.sp)[1]

## FORMAT DETECTOR LOCATIONS & NUMBER OF TRIALS PER DETECTOR IN ARRAYS/MATRICES
detector.xy <- coordinates(myDetectors$main.detector.sp)
n.trials <- as.vector(table(myDetectors$detector.sp$main.cell.id))
detector.dead.xy <- coordinates(myDetectors.dead$main.detector.sp)

## RETRIEVE DETECTION WINDOWS BOUNDARIES
lowerDetCoords <- detector.xy - 0.5 * myVars$DETECTORS$detResolution
upperDetCoords <- detector.xy + 0.5 * myVars$DETECTORS$detResolution
lowerDetCoords.dead <- detector.dead.xy - 0.5 * myVars$DETECTORS$detDeadResolution
upperDetCoords.dead <- detector.dead.xy + 0.5 * myVars$DETECTORS$detDeadResolution

## PLOT CHECK
if(myVars$plot.check){
   par(mfrow = c(1,2))
   ## PLOT NGS DETECTORS
   plot(myHabitat$buffered.habitat.poly, main = paste(n.detectors, "Detectors Alive"), col = rgb(0.16,0.67,0.16, alpha = 0.3))  
   plot(myStudyArea, add = TRUE, col = rgb(0.16,0.67,0.16,alpha = 0.5))
   plot(myDetectors$main.detector.sp, col = "red", pch = 16, cex = 0.1, add = TRUE)
   plot(COUNTRIES, add = TRUE)
   ## PLOT DEAD DETECTORS
   plot(myHabitat$buffered.habitat.poly, main = paste(n.detectors.dead, "Detectors Dead"), col = rgb(0.16,0.67,0.16, alpha = 0.3)) 
   plot(myStudyArea, add = T, col = rgb(0.16,0.67,0.16,alpha = 0.5))
   plot(myDetectors.dead$main.detector.sp, col = "red", pch = 16, cex = 0.1, add = TRUE)
   plot(COUNTRIES, add = TRUE)
}

### ====    3.2. GENERATE DETECTOR-LEVEL COVARIATES ====
### ====       3.2.1. EXTRACT COUNTRIES ====
dist <- gDistance(myDetectors$main.detector.sp, COUNTRIES, byid = T )
detCountries <- apply(dist,2, function(x) which.min(x))
detCountries <- as.numeric(as.factor(detCountries))

## PLOT CHECK 
if(myVars$plot.check){
   par(mfrow = c(1,2))
   myCol <- c("blue4", "yellow1")
   plot(GLOBALMAP, col = "gray80", main = "Countries")
   plot(myStudyArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
   plot(myBufferedArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
   points(myDetectors$main.detector.sp, col = myCol[detCountries], pch = 16, cex = 0.8)
   plot(COUNTRIES, add = TRUE)
}

### ====       3.2.2. EXTRACT COUNTIES ====
## ASSIGN COUNTIES TO DETECTORS
dist <- gDistance(myDetectors$main.detector.sp, COUNTIES_AGGREGATED , byid = TRUE)
detCounties <- apply(dist, 2, function(x) which.min(x))
COUNTIES_AGGREGATEDSubset <- COUNTIES_AGGREGATED[unique(detCounties),]
COUNTIES_AGGREGATEDSubset$idunique <- as.numeric(as.factor(unique(detCounties)))
detCounties <- as.numeric(as.factor(detCounties))
#detCounties <- COUNTIES_AGGREGATED$id[detCounties]

## PLOT CHECK 
if(myVars$plot.check){
   myCol <- terrain.colors(length(COUNTIES_AGGREGATED))
   plot(GLOBALMAP, col = "gray80", main = "Aggregated Counties")
   plot(myStudyArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
   plot(myBufferedArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
   points(myDetectors$main.detector.sp, col = myCol[detCounties], pch = 16, cex = 0.8)
   plot(COUNTIES_AGGREGATED, add = TRUE)
   text(COUNTIES_AGGREGATED, labels = COUNTIES_AGGREGATED$id, col = "black")  
}

### ====       3.2.3. EXTRACT GPS TRACKS LENGTHS ====
## INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detTracks <- matrix(NA, nrow = n.detectors, ncol = nYears)

## CREATE A LIST OF RASTERS
TRACKS.r <- list()
for(t in 1:nYears){
   ## EXTRACT THE NUMBER OF "TRACK-POINTS" PER RASTER CELL
   TRACKS.r[[t]] <- raster(myHabitat$buffered.habitat.poly, res = myVars$DETECTORS$detResolution)
   TRACKS.r[[t]][] <- 0
   tab1 <- table(cellFromXY(TRACKS.r[[t]], TRACKS_YEAR.sp[[t]]))
   TRACKS.r[[t]][as.numeric(names(tab1))] <- tab1
   detTracks[ ,t] <- raster::extract(TRACKS.r[[t]], myDetectors$main.detector.sp)
}#t


max <- max(unlist(lapply(TRACKS.r, function(x) max(x[]))))
cuts <- seq(0,max,length.out = 100)   #set breaks
col <- rev(terrain.colors(100))


# plot(habbdens, breaks=cuts, col = col,legend=FALSE, main=years[t]) #p
#    plot(myStudyArea.poly,add=T, border=grey(0.5))
#    points(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t],],
#           pch=16, cex=0.4, col=adjustcolor("black",alpha.f = 0.2))
#    plot(habbdens, legend.only=TRUE,breaks=cuts, col=col,
#         legend.width = 2,
#         axis.args=list(at=round(seq(0, max, length.out = 5),digits = 1),
#                        labels=round(seq(0, max, length.out = 5),digits = 1), 
#                        cex.axis=0.6),
#         legend.args=list(text='Density', side=4, font=2, line=2.5, cex=0.8))
#    
pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"Tracks.pdf",sep="")))

if(myVars$plot.check){
   # par(mfrow = c(2,2))#[CM]
   for(t in 1:nYears){
      plot( TRACKS.r[[t]],main=years[t], breaks=cuts, col = col,legend=FALSE)#[CM]
      plot(myHabitat$habitat.poly, main = years[t],add=T)
      # plot(myDetectors$main.detector.sp, cex = DoScale(detTracks[ ,t]), pch = 16, add = TRUE)
      plot( TRACKS.r[[t]], legend.only=TRUE,breaks=cuts, col=col,
                   legend.width = 2,
                   axis.args=list(at=round(seq(0, max, length.out = 5),digits = 1),
                                  labels=round(seq(0, max, length.out = 5),digits = 1),
                                  cex.axis=0.6),
                   legend.args=list(text='', side=4, font=2, line=2.5, cex=0.8))
      
      points(myFilteredData.sp$alive[myFilteredData.sp$alive$Year==years[t], ], col="red", pch=16, cex=0.8)
      
      }#t
}
dev.off()
### ====       3.2.4. EXTRACT DISTANCES TO ROADS ====
## AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- aggregate( DistAllRoads,
                           fact = myVars$DETECTORS$detResolution/res(DistAllRoads),
                           fun = mean)

## EXTRACT ROAD DISTANCE FOR EACH DETECTOR
detRoads <- raster::extract(DistAllRoads, myDetectors$main.detector.sp)

## if NA returns the average value of the cells within 20000m 
isna <- which(is.na(detRoads))
tmp <- raster::extract(DistAllRoads, myDetectors$main.detector.sp[isna,], buffer = 15000, fun = mean, na.rm = T)
detRoads[isna] <- tmp

if(myVars$plot.check){
   par(mfrow = c(1,1))
   plot(GLOBALMAP, col = "gray80", main = "Distance to roads")
   plot(myStudyArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
   plot(myBufferedArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
   plot(DistAllRoads,add=T)
   plot(myDetectors$main.detector.sp, cex=DoScale(detRoads), pch = 16, add = T)
}

### ====       3.2.5. EXTRACT DAYS OF SNOW ====
## EXTRACT SNOW 
detSnow <- matrix(0, nrow = dim(myDetectors$main.detector.sp)[1], ncol = nYears)
det.sptransf <- spTransform(myDetectors$main.detector.sp, CRS(proj4string(SNOW)))
detSnow[ ,1:nYears] <- raster::extract(SNOW, det.sptransf)

## if NA returns the average value of the cells within 20000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract(SNOW, det.sptransf[isna, ], buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:nYears] <- tmp
if(myVars$plot.check){
   plot(myDetectors$main.detector.sp,cex=DoScale(detSnow[,6],l = 0,u = 0.5),pch=16)
}
# still some NA... Increase buffer again 
isna <- which(is.na(detSnow),arr.ind = T)
isna <- unique(isna[,1])
tmp.list <- raster::extract(SNOW, det.sptransf[isna,], buffer = 35000)
detSnow[isna,1:nYears] <- unlist(lapply(tmp.list, function(x) colMeans(x, na.rm=T)))

### ====       3.2.6. SCALE AND ROUND DETECTOR-LEVEL COVARIATES ====
detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

detCovs <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2],3))
detCovs[,,1] <- detTracks
detCovs[,,2] <- matrix(detRoads,length(detRoads),nYears)
detCovs[,,3] <- detSnow

# CHECK IF CONTAINS NAs
if(any(is.na(detCovs))){print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")}



## PLOT CHECK
dev.off()

pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_detections over space and time.pdf",sep="")))
if(myVars$plot.check){
   for(t in 1:nYears){
      ## NGS DETECTIONS TOTAL
      tempTotal <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ]
      NGS_TabTotal <- table(tempTotal$Country)
      ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x) sum(x>0))
      ## ALIVE DETECTIONS INSIDE STUDY AREA/SAMPLING PERIOD
      tempIn <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ]
      NGS_TabIn <- table(tempIn$Country)
      ID_TabIn <- apply(table(tempIn$Id, tempIn$Country), 2, function(x) sum(x>0))
      ## PLOT NGS SAMPLES
      plot(GLOBALMAP, col="gray80")
      plot(myStudyArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
      plot(myBufferedArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
      # points(tempTotal, pch = 21, bg = "darkred")
      points(tempIn, pch = 21, bg = "blue")
      ## ADD NUMBER OF NGS samples and IDs per COUNTRY
      graphics::text(x = 100000, y = 7200000, labels = paste(NGS_TabTotal[names(NGS_TabTotal)=="N"],"NGS"), cex = 1.1, col = "firebrick3", font = 2)
      graphics::text(x = 100000, y = 7270000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
      graphics::text(x = 820000, y = 6780000, labels = paste(NGS_TabTotal[names(NGS_TabTotal)=="S"],"NGS"), cex = 1.1, col = "navyblue", font = 2)
      graphics::text(x = 820000, y = 6850000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
      ## ADD OVERALL NUMBERS
      mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
      mtext(text = paste(sum(NGS_TabIn), "NGS/", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
      mtext(text = paste(sum(NGS_TabTotal)-sum(NGS_TabIn), "NGS/", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"), side = 3, line = -1)
   }#t
}
dev.off()
### ====       1.2.2.deadData ====

## PLOT CHECK
if(myVars$plot.check){
   par(mfrow = c(1,3))
   for(t in 1:nYears){
      ## DEAD RECOVERIES TOTAL
      tempTotal <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
      NGS_TabTotal <- table(tempTotal$Country)
      ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x) sum(x>0))
      ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
      tempIn <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
      NGS_TabIn <- table(tempIn$Country)
      ID_TabIn <- apply(table(tempIn$Id, tempIn$Country), 2, function(x) sum(x>0))
      ## PLOT NGS SAMPLES
      plot(GLOBALMAP, col="gray80")
      plot(myStudyArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
      plot(myBufferedArea, col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
      # points(tempTotal, pch = 21, bg = "darkred")
      points(tempIn, pch = 21, bg = "blue")
      ## ADD NUMBER OF NGS samples and IDs per COUNTRY
      graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
      graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
      ## ADD OVERALL NUMBERS
      mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
      mtext(text = paste(sum(NGS_TabIn), "Dead Recoveries /", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
      # mtext(text = paste(sum(NGS_TabTotal), "Recoveries /", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"), side = 3, line = -1)
   }#t
}


## PLOT TREND DETECTIONS AND DEAD RECOVERIES OVER TIME AND SPACE #[CM]
#DETECTIONS
pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_TRENDDetections.pdf",sep="")))

temp <- unique(myFilteredData.sp$alive[,c("Year","Country","DNAID")]@data)
tab_Country.Year <- table(temp$Year, temp$Country)
country.colors <- c("goldenrod1","goldenrod3")

par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Detections", xlab="Years")
lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("bottomright",c("N","S"), fill=country.colors)

#ID DETECTED
temp <- unique(myFilteredData.sp$alive[,c("Year","Country","Id")]@data)
tab_Country.Year <- table(temp$Year, temp$Country)
country.colors <- c("goldenrod1","goldenrod3")

par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Id detected", xlab="Years")
lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("bottomright",c("N","S"), fill=country.colors)

## deadrecovery #[CM]
temp <- unique(myFilteredData.sp$dead.recovery[,c("Year","Country","Id")]@data)
tab_Country.Year <- table(temp$Year, temp$Country)
country.colors <- c("goldenrod1","goldenrod3")

par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)),
     ylab="N Id Dead recovered", xlab="Years")
lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("topright",c("N","S"), fill=country.colors)
dev.off()






### ==== 4.GENERATE y DETECTION ARRAYS ====
### ====    4.1.GENERATE NGS DETECTIONS : y.alive[i,j,t] ====
# allIDs <- unique(c(as.character(aliveData$Id), as.character(deadData$Id)))
# 
# ## CREATE NGS DETECTION ARRAYS
# source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/PD/Source/MakeY.R")
# source("C:/My_documents/RovQuant/Temp/PD/Source/MakeY.R")
# temp <- MakeY( samplesCoords = aliveData,
#                detectorCoords = cbind( coordinates(myDetectors$main.detector.sp),
#                                        myDetectors$main.detector.sp$main.cell.id),
#                subDetectorCoords = cbind( coordinates(myDetectors$detector.sp),
#                                           myDetectors$detector.sp$main.cell.id),
#                IDs = allIDs,
#                Years = NULL)
# 
# y.det.ALIVE <- temp$y
# numDets.ALIVE <- temp$numDets
# y.ar.ALIVE <- temp$y.binomial
# dimnames(y.ar.ALIVE) <- dimnames(temp$y.binomial)

## ASSIGN DETECTORS 
myData.alive <- AssignDetectors_v3( myData = myFilteredData.sp$alive
                                    ,                
                                    myDetectors = myDetectors$main.detector.sp
                                    ,
                                    mysubDetectors = myDetectors$detector.sp
                                    ,
                                    radius = myVars$DETECTORS$detResolution)

myData.dead <- AssignDetectors_v3( myData = myFilteredData.sp$dead.recovery
                                   ,
                                   myDetectors = myDetectors.dead$main.detector.sp
                                   ,
                                   radius = myVars$DETECTORS$detResolution)


##### EXPORT NGS DATA 
if(myVars$DATA$sex=="Hann"){
   assign("myFilteredData.spM", myFilteredData.sp)
   save(myFilteredData.spM, myFullData.spM,
        file = file.path(myVars$WD, myVars$modelName,
                         paste(myVars$modelName, "_NGSData", ".RData", sep = "")) )
   
   }else{
   assign("myFilteredData.spF", myFilteredData.sp)
      save(myFilteredData.spF, myFullData.spF,
           file = file.path(myVars$WD, myVars$modelName,
                            paste(myVars$modelName, "_NGSData", ".RData", sep = "")) )
      
}


## MAKE Y 
y.ar <- MakeY( myData = myData.alive$myData.sp
               ,
               myDetectors = myDetectors$main.detector.sp
               ,
               method = "Binomial"
               ,
               myData2 = myData.dead
               ,
               myDetectors2 = myDetectors.dead$main.detector.sp
               ,
               returnIdvector = TRUE
)


##PROJECT THE DEATH TO THE NEXT OCCASIONS 
y.ar.ALIVE <- y.ar$y.ar
dimnames(y.ar.ALIVE) <- dimnames(y.ar$y.ar)
dim(y.ar$y.ar)
table(myFilteredData.sp$alive$Year)
table(myData.alive$myData.sp$Year)
table(myData.dead$Year)


## PROJECT THE DEATH TO THE NEXT OCCASION.
y.ar.DEADProjected <- y.ar$y.ar2 
y.ar.DEADProjected[] <- 0
for(t in 2:nYears){y.ar.DEADProjected[,,t] <- y.ar$y.ar2[,,t-1]}

y.ar.DEAD <- apply(y.ar$y.ar2, c(1,3), function(x){if(sum(x)>0){which(x>0)}else{0}})
y.ar.DEAD <- cbind(rep(0, dim(y.ar.DEAD)[1]), y.ar.DEAD)
y.ar.DEAD <- y.ar.DEAD[ ,1:nYears]
dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])
dim(y.ar.DEAD)
y.ar.DEAD[y.ar.DEAD>0] <- 1
dim(y.ar$y.ar2)

## CHECK NGS DETECTIONS 
# if(myVars$plot.check){
#    par(mfrow = c(2,2))
#    ## sample some detected IDs randomly
#    detectedIDs <- which(apply(y.ar.ALIVE, 1, function(x)any(x>0)))
#    indexIDs <- sample(detectedIDs, 12)
#    namesIDs <- names(indexIDs)
#    for(i in 1:length(indexIDs)){
#       plot(myDetectors$main.detector.sp, pch = 19, cex = 0.5, col = "gray80", main = namesIDs[i])
#       plot(myHabitat$habitat.poly, add = TRUE)
#       plot(myStudyArea, add = TRUE)
#       ## plot activated detectors
#       detIndices <- unlist(apply(y.ar.ALIVE[indexIDs[i], , ],2,function(x)which(x > 0)))
#       points(detector.xy[detIndices,1], detector.xy[detIndices,2], pch = 19, cex = 1, col ="darkblue")
#       ## plot actual detection locations
#       plot(aliveData[aliveData$Id == namesIDs[i], ], add = T, col ="red")
#    }#i
# }
# 
# ### ====    4.2.GENERATE DEAD RECOVERIES : y.dead[i,t] ====
# ## CREATE DEAD RECOVERY DETECTION ARRAYS
# temp <- MakeY( samplesCoords = deadData,
#                detectorCoords = detector.dead.xy,
#                subDetectorCoords = NULL,
#                IDs = allIDs,
#                Years = NULL)
# 
# ## PROJECT THE DEATH TO THE NEXT OCCASION
# y.det.DEAD <- abind(array(0, c(dim(temp$y)[1],dim(temp$y)[2],dim(temp$y)[3],1)),
#                     array(temp$y[ , , ,-nYears], c(dim(temp$y)[1],dim(temp$y)[2],dim(temp$y)[3],dim(temp$y)[4]-1)),
#                     along = 4)
# dimnames(y.det.DEAD) <- dimnames(temp$y)
# numDets.DEAD <- cbind(rep(0,dim(temp$numDets)[1]), temp$numDets[ ,-nYears])
# y.ar.DEAD <- abind(array(0, c(dim(temp$y.binary)[1],dim(temp$y.binary)[2],1)), temp$y.binary[ , ,-nYears], along = 3)
# dimnames(y.ar.DEAD) <- dimnames(temp$y.binary)
# y.mx.DEAD <- apply(y.ar.DEAD, c(1,3), function(x){if(sum(x)>0){which(x>0)}else{0}})
# dimnames(y.mx.DEAD) <- list(dimnames(temp$y.binary)[[1]], dimnames(temp$y.binary)[[3]])
# 
# ## CHECK DEAD RECOVERIES
# if(myVars$plot.check){
#    par(mfrow = c(2,2))
#    recoveredIDs <- which(apply(y.mx.DEAD, 1, function(x)any(x>0)))
#    indexIDs <- sample(recoveredIDs, 12)
#    namesIDs <- names(indexIDs)
#    for(i in 1:length(indexIDs)){
#       plot(myDetectors.dead$main.detector.sp, pch = 19, cex = 0.5, col ="gray80")
#       plot(myHabitat$habitat.poly, add = TRUE)
#       plot(myStudyArea, add = TRUE)
#       points(detector.dead.xy[y.mx.DEAD[indexIDs[i],y.mx.DEAD[indexIDs[i], ] > 0],1],
#              detector.dead.xy[y.mx.DEAD[indexIDs[i],y.mx.DEAD[indexIDs[i], ] > 0],2],
#              pch = 19, cex = 1, col ="darkblue")
#       plot(deadData[deadData$Id == namesIDs[i], ], add = T, col ="red")
#    }#i
# }

## CONVERT TO BINARY DETECTIONS (recovered or not)
# y.mx.DEAD.binary <- y.mx.DEAD
# y.mx.DEAD.binary[y.mx.DEAD.binary > 0] <- 1

### ====    4.3.CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ====
distances <- list()
for(t in 1:nYears){
   print(paste("------ ", t ," -------", sep = "" ))
   distances[[t]] <- CheckDistanceDetectionsV2( y = y.ar.ALIVE[,,t], 
                                                detector.xy = coordinates(myDetectors$main.detector.sp), 
                                                max.distance = myVars$DETECTIONS$maxDetDist,
                                                method = "pairwise",
                                                plot.check = T)
      
      
   # PLOT INDIVIDUALS THAT DO HAVE DETECTIONS FURTHER AWAY THAN THRESHOLD DISTANCE
   if(myVars$plot.check){
      par(mfrow = c(1,1))
      if(sum(distances[[t]]$y.flagged) > 0){
         affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
         for(i in affected.ids){
            plot(myStudyArea, main = paste("t: ",t,"     i: ", i, sep = ""))
            scalebar(2*myVars$DETECTIONS$maxDetDist, xy = c(800000,6700000), type = "bar", divs = 2, below = "km",
                     label = c(0, myVars$DETECTIONS$maxDetDist/1000, myVars$DETECTIONS$maxDetDist/500), cex = 0.8, adj = c(0.5,-0.9))
            plot(COUNTRIES, add = T)
            plot(myDetectors$main.detector.sp, add = T, col = grey(0.8), cex = 0.3, pch = 19)
            
            tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Id == dimnames(y.ar.ALIVE)[[1]][i] & myFilteredData.sp$alive$Year == years[t], ]
            tmp <- tmp[order(tmp$Date), ]
            tmp.xy <- coordinates(tmp)
            n.det <- nrow(tmp.xy)
            
            points(tmp, col = "pink", pch = 16, cex = 1)
            arrows(x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
                   x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2],
                   length = 0.1, lwd = 1)
            points(myDetectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0), ], pch = 16, col = "red")
            
            tmp2 <- myDetectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0 & distances[[t]]$y.flagged[i,] == 1), ]
            plot(tmp2, add = T, col = "blue", pch = 13, cex = 1.5, lwd = 1)
         }#i
      }#if
   }#if plot.check
   ## REMOVE DETECTIONS THAT ARE FURTHER THAN  THE THRESHOLD
   y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
}#t

### ====    4.4.GENERATE INDIVIDUAL-LEVEL COVARIATES ====
### ====       4.4.1.TRAP-RESPONSE ====
## Make matrix of previous capture indicator
# already.detected <- MakeTrapResponseCov(data = mycleanData, IDs = allIDs)
# ## Subset to focal years
# already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar.ALIVE)[[3]]]
# ## Subset to focal individuals
# already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar.ALIVE)[[1]], ]
## Make matrix of previous capture indicator
already.detected <- MakeTrapResponseCov(myFullData.sp$alive, myFullData.sp$dead.recovery)

## Subset to focal years
already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar.ALIVE)[[3]]]

## Subset to focal individuals
already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar.ALIVE)[[1]], ]



## Plot an image of the matrix
if(myVars$plot.check){
   par(mfrow = c(1,1))
   barplot(colSums(apply(y.ar.ALIVE, c(1,3),function(x)any(x>0))))
   barplot(colSums(already.detected), add = TRUE, col = "gray40")
   legend(x = 0, y = 250, legend = c("newly Det", "already Det"),
          fill = c("gray80", "gray40"))
}

### ====       4.4.2.TELEPORTATION COVARIATE ====
# source("C:/My_documents/RovQuant/Temp/PD/Source/CheckDistanceACS.R")
# source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/PD/Source/CheckDistanceACS.R")
# 
# distancesACs <- CheckDistanceACS( y = y.ar.ALIVE,
#                                   detector.sp = myDetectors$main.detector.sp,
#                                   myStudyArea.sp = myStudyArea,
#                                   maxDist = myVars$DETECTIONS$maxMoveDist,
#                                   plot.check = myVars$plot.check)
# 
# long.dispersal.id <- which(distancesACs$distances.not.consecutive > myVars$DETECTIONS$maxMoveDist, arr.ind = TRUE)
# 
# dispersalToggle <- matrix(0, nrow = dim(y.ar.ALIVE)[1], ncol = nYears)
# for(i in 1:nrow(long.dispersal.id)){
#    dispersalToggle[long.dispersal.id[i,2],long.dispersal.id[i,2]] <- 1
# }#i   


### ====    5.4. AGE ====#[CM]
min.age <- age <- precapture <- matrix(NA, dim(y.ar.ALIVE)[1], dim(y.ar.ALIVE)[3], dimnames = list(y.ar$Id.vector,years))

temp <- apply(y.ar.ALIVE, c(1,3), sum)
year.first.capture <- apply(temp, 1, function(x)min(years[which(x>0)]))
year.first.capture[is.infinite(year.first.capture)] <- NA
names(year.first.capture) <- y.ar$Id.vector

for(i in y.ar$Id.vector){
   this.set <- myData.dead[myData.dead$Id == i, ]
   year.dead <- myData.dead$Death[myData.dead$Id == i]
   year.first.captured <- year.first.capture[i]
   precapture[i,] <- as.numeric(years < year.first.captured)
   if(all(is.na(precapture[i,])))precapture[i,] <- 1
   latest.recruitment.year <- min(year.dead,year.first.captured, na.rm = TRUE) 
   
   try({
      min.age[i,] <- years-latest.recruitment.year
   },silent = TRUE)
   
   try({
      birth.year <- this.set@data$Death-this.set@data$min.age
      if(birth.year<latest.recruitment.year) min.age[i,] <- years-birth.year 
   },silent = TRUE)
   
   try({
      birth.year <- this.set@data$Death - this.set@data$age
      age[i,] <- years-birth.year
   }, silent = TRUE)
}
image(t(min.age))
image(t(age))


### ==== 6. CONVERT TO CACHED DETECTORS AND SPARSE MATRIX ====
### ====   6.1 RESCALE COORDINATES  ====
# HABITAT
ScaledLowerCoords <- UTMToGrid(data.sp = SpatialPoints(lowerHabCoords),
                               grid.sp = myHabitat$habitat.sp)$data.scaled.xy
ScaledUpperCoords <- UTMToGrid(data.sp = SpatialPoints(upperHabCoords),
                               grid.sp = myHabitat$habitat.sp)$data.scaled.xy
ScaledUpperCoords[ ,2] <- ScaledUpperCoords[ ,2]+1
ScaledLowerCoords[ ,2] <- ScaledLowerCoords[ ,2]-1

# DETECTORS
ScaledDetectors <- UTMToGrid(data.sp = myDetectors$main.detector.sp,
                             grid.sp = myHabitat$habitat.sp)$data.scaled.xy

# SXY
nimInits$sxy <- UTMToGrid(data.sxy = sxy.init,
                          grid.sp = myHabitat$habitat.sp)$data.scaled.xy
nimData$sxy <- UTMToGrid(data.sxy = nimData$sxy,
                         grid.sp = myHabitat$habitat.sp)$data.scaled.xy
# PLOT CHECK 
plot(ScaledUpperCoords[ ,2] ~ ScaledUpperCoords[ ,1])
points(ScaledLowerCoords[ ,2] ~ ScaledLowerCoords[ ,1], col = "red")
plot(ScaledUpperCoords[ ,2] ~ ScaledUpperCoords[ ,1], xlim = c(0,5), ylim = c(0,5))
points(ScaledLowerCoords[ ,2] ~ ScaledLowerCoords[ ,1], col = "red")
points(nimData$detector.xy[ ,2] ~ nimData$detector.xy[ ,1], col = "blue")

for(t in 1:nimConstants$n.years){
  points(nimInits$sxy[,2,t] ~ nimInits$sxy[,1,t], col = "green", pch = 16)
}

# ADD TO NIMDATA
nimData$detector.xy <- as.matrix(ScaledDetectors)          
nimData$lowerHabCoords <- as.matrix(ScaledLowerCoords)
nimData$upperHabCoords <- as.matrix(ScaledUpperCoords)

### ====   6.2 CREATE CACHED DETECTORS OBJECTS ====
#[CM] reduce multiplicator to 3 
maxDistReCalc <- 2.1*myVars$DETECTIONS$maxDetDist #+ sqrt(2*(myVars$DETECTIONS$resizeFactor*myVars$HABITAT$habResolution)^2)

DetectorIndexLESS <- GetDetectorIndexLESS( habitat.mx = myHabitat$habitat.mx,
                                           detectors.xy = nimData$detector.xy,
                                           maxDist = maxDistReCalc/res(myHabitat$habitat.r)[1],
                                           ResizeFactor = 2,
                                           plot.check = TRUE)
DetectorIndexLESS$nDetectorsLESS
# ADD TO NIMDATA
dim(DetectorIndexLESS$detectorIndex)
nimConstants$y.maxDet <- dim(DetectorIndexLESS$habitatID)[1]
nimConstants$x.maxDet <- dim(DetectorIndexLESS$habitatID)[2]
nimConstants$ResizeFactor <- DetectorIndexLESS$ResizeFactor
nimConstants$n.cellsSparse <- dim(DetectorIndexLESS$detectorIndex)[1]
nimConstants$maxNBDets <- DetectorIndexLESS$maxNBDets

nimData$detectorIndex <- DetectorIndexLESS$detectorIndex
nimData$nDetectorsLESS <- DetectorIndexLESS$nDetectorsLESS
nimData$habitatIDDet <- DetectorIndexLESS$habitatID

### ====   6.3 TRANSFORM Y TO SPARSE MATRICES  ====
SparseY <- GetSparseY(y.alive)

# ADD TO NIMDATA
nimData$y.alive <- SparseY$y 
nimData$yDets <- SparseY$yDets
nimData$nbDetections <- SparseY$nbDetections
nimConstants$nMaxDetectors <- SparseY$nMaxDetectors



### ==== 5. MAKE AUGMENTATION ====
## DATA ARRAYS
y.alive <- MakeAugmentation(y = y.ar.ALIVE, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.dead <- MakeAugmentation(y = y.ar.DEAD, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
# y.dead.binary <- MakeAugmentation(y = y.mx.DEAD.binary, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)

## INDIVIDUAL COVARIATES
already.detected <- MakeAugmentation(y = already.detected, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
# dispersalToggle <- MakeAugmentation(y = dispersalToggle, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
age <- MakeAugmentation(y = age, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
min.age <- MakeAugmentation(y = min.age, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
precapture <- MakeAugmentation(y = precapture, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)

## -----------------------------------------------------------------------------
## ------ III.MODEL SETTING AND RUNNING ------- 
### ==== 1. NIMBLE MODEL DEFINITION ====
modelCode <- nimbleCode({
   ##--------------------------------------------------------------------------------------------
   ##-----------------------------## 
   ##------ SPATIAL PROCESS ------##  
   ##-----------------------------##  
   # dispSigma ~ dgamma(0.01,0.01)
   # betaDens ~ dnorm(mean = 0, sd = 5)
   #[CM]
   dispSigma ~ dunif(0,10)
   betaDens  ~ dnorm(0.0,0.01)
   # [CM]
   # for(h in 1:numHabWindows){
   #    log(mu[h]) <- betaDens * denCounts[h]
   # }#h
   #mu[1:numHabWindows] <- exp(betaDens * denCounts[1:numHabWindows,1])
   
   
   habIntensity[1:numHabWindows] <- exp(betaDens * denCounts[1:numHabWindows,1])
   sumHabIntensity <- sum(habIntensity[1:numHabWindows])
   logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
   logSumHabIntensity <- log(sumHabIntensity)
   
   for(i in 1:n.individuals){
      sxy[i, 1:2, 1] ~ dbernppAC(
         lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
         upperCoords = upperHabCoords[1:numHabWindows, 1:2],
         logIntensities = logHabIntensity[1:numHabWindows],
         logSumIntensity = logSumHabIntensity,
         habitatGrid = habitatGrid[1:y.max,1:x.max],
         numGridRows =  y.max,
         numGridCols = x.max
         
      )
   }#i
   
  for(t in 2:n.years){
      for(i in 1:n.detected){
         
         
         sxy[i, 1:2, t] ~ dbernppACmovement_normal(
            lowerCoords = lowerHabCoords[1:numHabWindows, 1:2]
            ,
            upperCoords = upperHabCoords[1:numHabWindows, 1:2]
            ,
            s = sxy[i, 1:2, t - 1]
            ,
            sd = dispSigma
            ,
            baseIntensities = habIntensity[1:numHabWindows]
            ,
            habitatGrid =  habitatGrid[1:y.max,1:x.max]
            ,
            numGridRows = y.max
            ,
            numGridCols = x.max
            ,
            numWindows= numHabWindows
         )
         
         
      }#i
      
      for(i in n.detected1:n.individuals){
         # sxy[i,1:2,t] ~ dbinomPPSingle( lowerHabCoords[1:numHabWindows,1:2],
         #                                upperHabCoords[1:numHabWindows,1:2],
         #                                mu[1:numHabWindows],
         #                                1, numHabWindows)
         
         sxy[i, 1:2, t] ~ dbernppAC(
            lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
            upperCoords = upperHabCoords[1:numHabWindows, 1:2],
            logIntensities = logHabIntensity[1:numHabWindows],
            logSumIntensity = logSumHabIntensity,
            habitatGrid = habitatGrid[1:y.max,1:x.max],
            numGridRows =  y.max,
            numGridCols = x.max
            
         )
         
      }#i
   }#t
   
   ##--------------------------------------------------------------------------------------------
   ##-------------------------------## 
   ##----- DEMOGRAPHIC PROCESS -----## 
   ##-------------------------------##    
   omeg1[1:3] ~ ddirch(alpha[1:3])   
   
   for(t in 1:n.years1){
      # PRIORS 
      gamma[t] ~ dunif(0,1)
      w[t] ~ dunif(0,0.59)
      h[t] ~ dunif(0,0.39)
      phi[t] <- 1-h[t]-w[t]
      
      # "UNBORN"
      omega[1,1,t] <- 1-gamma[t]
      omega[1,2,t] <- gamma[t]
      omega[1,3,t] <- 0
      omega[1,4,t] <- 0
      # "NON-PAIRS"
      omega[2,1,t] <- 0
      omega[2,2,t] <- phi[t]
      omega[2,3,t] <- h[t]
      omega[2,4,t] <- w[t]
      # "NEWLY DEAD LEGAL HUNTING"[RB]
      omega[3,1,t] <- 0
      omega[3,2,t] <- 0
      omega[3,3,t] <- 0
      omega[3,4,t] <- 1
      # "NEWLY DEAD OTHER SOURCES AND DEAD"[RB]
      omega[4,1,t] <- 0
      omega[4,2,t] <- 0
      omega[4,3,t] <- 0
      omega[4,4,t] <- 1
   }#t
   
   for(i in 1:n.individuals){ 
      z[i,1] ~ dcat(omeg1[1:3]) 
      for(t in 1:n.years1){
         z[i,t+1] ~ dcat(omega[z[i,t],1:4,t]) 
      }#i 								
   }#t 
   
   ##---------------------------------------------------------------------------------------------   
   ##-----------------------------##
   ##----- DETECTION PROCESS -----## 
   ##-----------------------------##
   sigma ~ dunif(0,4)
   
   # betaCovs[1] <- betaTracks
   # betaCovs[2] <- betaRoads
   # betaCovs[3] <- betaSnow
   for(c in 1:n.covs){
      betaCovs[c] ~ dunif(-5,5)
   }
   
   betaResponse ~ dunif(-5,5)
   
   for(c in 1:n.countries){
      for(t in 1:n.years){
         p0[c,t] ~ dunif(0,1)
      }#t
   }#c  
   
   # for(c in 1:n.countries){
   #    for(t in 1:n.years){
   #       p0[c,1,t] ~ dunif(0,1)
   #       logit(p0[c,2,t]) <- logit(p0[c,1,t]) + betaResponse
   #    }#t
   # }#c     
   
   for(t in 1:n.years){
      for(i in 1:n.individuals){
         # y.alive[i,1:nMaxDetectors,t] ~ dbin_LESS_Cached_MultipleCov( sxy = sxy[i,1:2,t],
         #                                                              sigma = sigma,
         #                                                              nbDetections[i,t],
         #                                                              yDets = yDets[i,1:nMaxDetectors,t],
         #                                                              detector.xy =  detector.xy[1:n.detectors,1:2],
         #                                                              trials = trials[1:n.detectors],
         #                                                              detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
         #                                                              nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
         #                                                              ResizeFactor = ResizeFactor,
         #                                                              maxNBDets = maxNBDets,
         #                                                              habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
         #                                                              indicator = isAlive[i,t],
         #                                                              p0[1:n.countries,detResponse[i,t],t],
         #                                                              detCountries[1:n.detectors],
         #                                                              detCov = detCovs[1:n.detectors,t,1:n.covs],
         #                                                              betaCov = betaCovs[1:n.covs])
         y.alive[i,1:nMaxDetectors,t] ~ dbin_LESS_Cached_MultipleCovResponse(  sxy = sxy[i,1:2,t],
                                                                               sigma = sigma,
                                                                               nbDetections[i,t],
                                                                               yDets = yDets[i,1:nMaxDetectors,t],
                                                                               detector.xy =  detector.xy[1:n.detectors,1:2],
                                                                               trials = trials[1:n.detectors],
                                                                               detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
                                                                               nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
                                                                               ResizeFactor = ResizeFactor,
                                                                               maxNBDets = maxNBDets,
                                                                               habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                               indicator = isAlive[i,t],
                                                                               p0[1:n.countries,t],
                                                                               detCountries[1:n.detectors],
                                                                               detCov = detCovs[1:n.detectors,t,1:n.covs],
                                                                               betaCov = betaCovs[1:n.covs],
                                                                               BetaResponse = betaResponse,
                                                                               detResponse = detResponse[i,t])
         
         y.dead[i,t] ~ dbern(z[i,t] == 3) 
      }#i
   }#t
   
   ##---------------------------------------------------------------------------------------------										  
   ##----------------------------------------## 
   ##---------- DERIVED PARAMETERS ----------##
   ##----------------------------------------##
   for(i in 1:n.individuals){ 
      isAlive[i,1] <- (z[i,1] == 2) 
      for(t in 1:n.years1){
         isAlive[i,t+1] <- (z[i,t+1] == 2) 
      }
   }
   for(t in 1:n.years){
      N[t] <- sum(isAlive[1:n.individuals,t])
   }#t
})

### ==== 2. NIMBLE CONSTANTS ====
nimConstants <- list( n.individuals = dim(y.alive)[1],
                      n.detected = dim(y.ar.ALIVE)[1],
                      n.detected1 = dim(y.ar.ALIVE)[1]+1,
                      n.detectors = dim(y.alive)[2],  
                      n.years = dim(y.alive)[3], 
                      n.years1 = dim(y.alive)[3]-1, 
                      n.covs = dim(detCovs)[3],
                      numHabWindows = nHabCells,
                      n.countries = max(detCounties),
                      y.max = dim(habIDCells.mx)[1],
                      x.max = dim(habIDCells.mx)[2])#max(detCountries))

### ==== 3. NIMBLE INITS ====
### ====    3.1.GENERATE INITIAL z ====
   y.init <- apply(y.alive,c(1,3), function(x)as.numeric(any(x > 0)))
   # y.init <- y.init + 2*y.dead.binary + 1
   y.init <- y.init + 2*y.dead + 1
   
   STATE <- matrix(c(1,1,0,0,
                     0,1,1,1,
                     0,0,0,1,
                     0,0,0,1), nrow = 4, byrow = TRUE)
   OBS <- matrix(c(1,0,0,
                   1,1,0,
                   0,0,1,
                   1,0,0), nrow = 4, byrow = TRUE)

   z.init <- MakeZ( y = y.init,
                    STATE = STATE,
                    OBSERVATION = OBS,
                    f.state = c(1,2),
                    z.init.NA = TRUE)
   
   # ## Reconstruct monthly z based on ALL detections and dead recoveries
   # zMonths <- MakeZfromScratch( data.alive = myFullData.sp$alive,
   #                              data.dead = myFullData.sp$dead.recovery,
   #                              samplingMonths = unlist(myVars$DATA$samplingMonths))
   # table(zMonths)
   # 
   # 
   # 
   # ## Subset to focal years
   # zMonths <- zMonths[ , ,dimnames(zMonths)[[3]] %in% dimnames(y.alive)[[3]]]
   # 
   # 
   # image(t(zMonths[,,1]))
   # 
   # ## Subset to focal individuals
   # zMonths <- zMonths[dimnames(zMonths)[[1]] %in% dimnames(y.alive)[[1]], , ]
   # 
   # ## Augment z1 
   # zMonths <- MakeAugmentation(y = zMonths,
   #                             aug.factor = myVars$DETECTIONS$aug.factor,
   #                             replace.value = NA)
   # 
   # ## Compress back to yearly z
   # zYears <- apply(zMonths, c(1,3), function(x){
   #    if(any(x[1:length(unlist(myVars$DATA$samplingMonths))] == 1, na.rm = T)){
   #       2
   #    }else{
   #       if(any(x[1:length(unlist(myVars$DATA$samplingMonths))] >= 2, na.rm = T)){
   #          4
   #       }else{
   #          0
   #       }}})
   # 
   # zYears[60,]
   # 
   # table(zYears[])
   # ## Combine with the Social State information
   # indSocialState<-1
   # zYears <- zYears + indSocialState
   # zYears[zYears <= 2] <- NA
   # zYears[zYears == 3] <- 2
   # zYears[zYears == 4] <- 3
   # zYears[zYears >= 5] <- 4
   # table(zYears)
   # 
   # 
   # 
   # 
   # ## INCORPORATE AGE
   # 
   # z <- zYears
   # z.age <- z
   # 
   ### ====    3.2.GENERATE INITIAL sxy ====
   #[CM]
   # sxy.init <- MakeSxyInits(y = y.alive,
   #                          y.dead = y.dead,
   #                          habitatPoly = myHabitat$buffered.habitat.poly,
   #                          detCoords = detector.xy,
   #                          detNum = NULL,
   #                          detCoords.dead = detector.dead.xy,
   #                          maxDetDist = myVars$DETECTIONS$maxDetDist,
   #                          maxMoveDist = NULL,
   #                          plot.check = F)
   #remove an annoying detection
   # i=443
   # t=3
   # y.alive[i,which(y.alive[i,,t]>0)[1],t] <- 0
                            

   sxy.init <- MakeInitsXY( y = y.alive,
                            detector.xy = detector.xy,
                            habitat.r = myHabitat$habitat.r,
                            ydead = y.dead,
                            detector.xyDead = detector.dead.xy,
                            dist.move = 10000)
   
   ## ADD SXY DATA 
   habitat.poly <- aggregate(rasterToPolygons(myHabitat$habitat.r,fun = function(x){x==1}))
   proj4string(habitat.poly) <- CRS(proj4string(myHabitat$habitat.sp))
   z.age <- z.init$z.reconstruct.mx
   sxy.data <- sxy.init
   sxy.data[] <- NA
   for(i in 1:length(y.ar$Id.vector)){
      for(t in 1:(dim(sxy.data)[3]-1)){
         if(sum(z.age[i,t+1] %in% c(3,4))>0 ){
            temp <- myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Id == y.ar$Id.vector[i] & myFullData.sp$dead.recovery$Year == years[t], ]
            if(length(temp) > 0){
               if(!is.na(over(temp,habitat.poly))){
                  if(raster::extract(myHabitat$habitat.r, temp)==0){
                     #id could be dead in the spatial extent but can be outside the habitat because in a cell <49%habitat
                     buff <- gBuffer(temp,width = myHabitat$resolution*2)
                     inter <- raster::intersect(buff,habitat.poly)
                     sxy.data[i, ,t+1]  <- coordinates(spsample(x = inter,n = 1,type="random"))
                  }else{sxy.data[i, ,t+1] <- coordinates(temp)}
               }
            }
         }
      }
   }
   sxy.init[!is.na(sxy.data)] <- NA
   
   which(!is.na(sxy.data),arr.ind = T)
   
   ### ====    3.4. LIST NIMBLE INITS ====
   nimInits <- list( "sxy" = sxy.init/1000,
                     "dispSigma" = runif(1,0,10),
                     "z" = z.init$z.init.mx,
                     "omeg1" = c(0.5,0.25,0.25),
                     "gamma" = runif(dim(y.alive)[3]-1,0,1),
                     "p0" = array(runif(18,0,0.2), c(nimConstants$n.countries,dim(y.alive)[3])),
                     "sigma" = runif(1,0,4),
                     "betaDens" = runif(1,-0.1,0.1),#[CM]#0,
                     "betaCovs" = runif(dim(detCovs)[3],-0.1,0.1),#[CM]rep(0,dim(detCovs)[3]),
                     "betaResponse" = runif(1,-0.1,0.1),#[CM]#0,
                     "h" = runif(dim(y.alive)[3]-1,0.1,0.3), #  runif(dim(y.alive)[3]-1,0.2,0.4),
                     "w" = runif(dim(y.alive)[3]-1,0.1,0.3)) #,runif(dim(y.alive)[3]-1,0.2,0.4))
   
 
   
   ### ==== 4. NIMBLE DATA ====
   nimData <- list( z = z.init$z.reconstruct.mx,   
                    y.alive = y.alive,
                    y.dead = y.dead,#y.dead.binary,
                    lowerHabCoords = lowerHabCoords/1000, 
                    upperHabCoords = upperHabCoords/1000, 
                    detCountries = detCounties,#detCountries,
                    detCovs = detCovs,
                    detResponse = already.detected ,#+ 1,
                    denCounts = denCounts,
                    trials = n.trials,
                    alpha = rep(1,3),
                    detector.xy = detector.xy/1000,
                    sxy = sxy.data,
                    habitatGrid = habIDCells.mx )
   
   ### ==== 5. NIMBLE PARAMETERS ====
   nimParams <- c("N", 
                  "dispSigma",
                  "omeg1", "gamma", "phi", "h", "w",
                  "p0", "sigma", "betaDens", "betaCovs","betaResponse",
                  "z", "sxy")
   
   ### ==== 7. SAVE NIMBLE INPUT ====
   save(nimData,
        nimConstants,
        nimParams,
        modelCode,
        nimInits,
        file = file.path(myVars$WD, myVars$modelName,
                         paste(myVars$modelName, c, ".RData", sep = "")))


apply(nimData$z, 2, function(x) sum(x==2,na.rm = T))



### ==== 8. NIMBLE RUN ====
load(file.path(myVars$WD, myVars$modelName,
               paste(myVars$modelName, c, ".RData", sep = "")))

model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      inits = nimInits,
                      data = nimData,
                      check = FALSE,
                      calculate = FALSE) 

conf <- configureMCMC(model, monitors = nimParams, thin = 1)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = model, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc



### ====    8.4. RUN NIMBLE MCMC IN SUCCESSIVE BITES ====
## SET NUMBER OF BITES AND NUMBER OF ITERATIONS PER BITE
bite.size <- 100 
bite.number <- 10

## LOOP OVER NUMBER OF BITES
for(nb in 1:bite.number){
   print(nb)
   if(nb == 1){
      ## run initial MCMC
      MCMCRuntime <- system.time(Cmcmc$run(bite.size))
   } else {      
      ## run subsequent MCMCs
      MCMCRuntime <- system.time(Cmcmc$run(bite.size, reset = FALSE))
   }
   
   ## STORE BITE OUTPUT IN A MATRIX
   mcmcSamples <- as.matrix(Cmcmc$mvSamples)
   CumulRunTime <- proc.time() - ptm
   
   ## EXPORT NIMBLE OUTPUT 
   outname <- file.path(path.OUT, paste("NimbleBite", nb, "_FOR", set, sep = ""))
   save(CumulRuntime, MCMCRuntime, mcmcSamples, file = outname)
   
   ## FREE UP MEMORY SPACE 
   rm("mcmcSamples") 
   Cmcmc$mvSamples$resize(0) ## reduce the internal mvSamples object to 0 rows,
   gc() ## run R's garbage collector
}#nb

## -----------------------------------------------------------------------------
## ------ IV.PROCESS RESULTS ------
### ==== 1. LOAD & PROCESS NIMBLE OUTPUTS ====
# List the directories containing bite outputs
myVars$modelName <- myVars$modelName#"23.J_Fa"
outDirectories <- list.files(file.path(myVars$WD, myVars$modelName))[grep("NimbleOut", list.files(file.path(myVars$WD, myVars$modelName)))]
path.list <- file.path(myVars$WD, myVars$modelName, outDirectories)

# Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x){
   files <- list.files(x)
   files <- files[grep(".RData", files)]
   length(files)
}))
minBites <- min(numBites)


#Niter to Remove (burn-in) #CM
NSkipBites <- 18
nthin <- 1
nimOutput <- RUNTIME <- list()
gc()
for(p in 1:length(path.list)){
   print(path.list[p])
   outfiles <- list.files(path.list[p])
   out <- runtime <- list()#[CM]
   for(x in NSkipBites:minBites){
      print(x)
      load(file.path(path.list[p], paste("bite_", x, ".RData", sep = "")))
      runtime[[x]] <- RunTime[3] 
      params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
      parmIndex <- which(! params.simple %in% c("sxy","z"))
      nthins <- seq(1,dim(this.sample)[1], by=nthin)
      out[[x]] <- this.sample[nthins,parmIndex]#[ ,parmIndex] 
   }#x
   RUNTIME[[p]] <- unlist(runtime)#[CM]
   out.mx <- do.call(rbind, out)
   nimOutput[[p]] <- as.mcmc(out.mx)
}#p


lapply(RUNTIME, function(x) x/3600)#[CM]
unlist(lapply(RUNTIME, function(x) x/3600))#[CM]
TIME <- lapply(RUNTIME, function(x) x/3600)

max <- unlist(lapply(RUNTIME, function(x) x/3600))
at=c(1:length(TIME[[1]]))
plot(TIME[[1]]~at, pch=16, col=adjustcolor("red", alpha.f = 0.5), xlim=c(0,length(TIME[[1]])+2),
     ylim=c(0,10), ylab="time hours", xlab="bite number")
for(i in 2:length(TIME)){
   points(TIME[[i]] ~ at, )
}
nimOutput <- as.mcmc.list(nimOutput)
myResults <- ProcessCodaOutput(nimOutput,params.omit = c("sxy","z"))



for(i in 1: length(nimOutput)){

   nimOutput[[i]][,"sigma"] <-  nimOutput[[i]][,"sigma"]*myHabitat$resolution
   nimOutput[[i]][,"dispSigma"] <-  nimOutput[[i]][,"dispSigma"]*myHabitat$resolution
}

{#doall[CM]



   
   
   
   ### ==== 2. PLOT PARAMETERS ESTIMATES ====
   pdf(file = file.path( myVars$WD,
                         myVars$modelName,paste(myVars$modelName,".pdf",sep="")
   ))
   
   ### ====    2.1.PLOT DETECTIONS ====
   myLayOut.mx <- cbind(c(2,1), c(1,1))
   myLayOut <- layout(myLayOut.mx, width = c(1,2), heights = c(1,2))
   #layout.show(myLayOut)
   
   ## PLOT STUDY AREA
   par(mar = c(0,0,0,0))
   plot(myHabitat$buffered.habitat.poly,  col = rgb(34/250, 139/250, 34/250, alpha = 0))
   plot(GLOBALMAP, col = "gray80", add = TRUE)
   plot(myHabitat$habitat.poly, col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
   plot(myHabitat$buffered.habitat.poly, col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
   
   ## PLOT DETECTIONS
   points(myFilteredData.sp$alive, pch = 3, col = "darkred", cex = 0.5)
   points(myFilteredData.sp$dead.recovery, pch = 3, col = "darkblue", cex = 0.5)
   
   ## ADD NUMBER OF DETECTIONS
   graphics::text(x = 190000, y = 7820000, cex = 1.5,
                  labels = paste(dim(myFilteredData.sp$alive)[1], "NGS samples"))
   
   ## ADD NUMBER OF INDIVIDUALS DETECTED
   graphics::text(x = 190000, y = 7880000, cex = 1.5,
                  labels = paste(length(unique(myFilteredData.sp$alive$Id)), "Individuals"))
   
   # ## PLOT GLOBAL MAP
   # par(mar = c(0,0,0,0))
   # plot(GLOBALMAP, col = "gray80")
   # plot(myStudyArea, col = "red", add = TRUE)
   
   ### ====    2.1.N ==== 
   par(mfrow = c(1,1), mar = c(5,5,5,5))
   plot(10, xlim = c(0, nYears+1), ylim = c(200,1200), type ="n", xaxt="n", xlab = "Years", ylab = "N")
   axis(1, c(1:nYears),labels = years)
   for(t in 1:nYears){
      plot.violins(list(myResults$sims.list$N[,t]),
                   x = t,
                   at = t,
                   violin.width = 0.3,
                   col = "firebrick3",
                   add = T,
                   alpha = 0.3,
                   border.col = "firebrick3")
   }#t
   
   params <- dimnames(nimOutput[[1]])[[2]][grep("N",dimnames(nimOutput[[1]])[[2]])]
   for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
   }
   
   ### ====    2.2.h ==== 
   par(mfrow = c(1,1))
   plot(10, xlim = c(0, nYears), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "h")
   axis(1, c(1:nYears),labels = years)
   myCol <- "firebrick3"
   for(t in 1:(nYears-1)){
      plot.violins(list(myResults$sims.list$h[ ,t]),
                   x = t,
                   at = t,
                   violin.width = 0.3,
                   col = myCol,
                   add = T,
                   alpha = 0.3,
                   border.col = myCol)
   }#t
   
   params <- dimnames(nimOutput[[1]])[[2]][grep("h\\[",dimnames(nimOutput[[1]])[[2]])]
   for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
   }
   
   ### ====    2.3.gamma ==== 
   par(mfrow = c(1,1))
   plot(10, xlim = c(0, nYears), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "gamma")
   axis(1, c(1:nYears),labels = years)
   myCol <- "firebrick3"
   for(t in 1:(nYears-1)){
      plot.violins(list(myResults$sims.list$gamma[ ,t]),
                   x = t,
                   at = t,
                   violin.width = 0.3,
                   col = myCol,
                   add = T,
                   alpha = 0.3,
                   border.col = myCol)
   }#t
   
   params <- dimnames(nimOutput[[1]])[[2]][grep("gamma",dimnames(nimOutput[[1]])[[2]])]
   for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
   }
   
   ### ====    2.4.w ==== 
   par(mfrow = c(1,1))
   plot(10, xlim = c(0, nYears+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "w")
   axis(1, c(1:nYears),labels = years)
   for(t in 1:(nYears-1)){
      plot.violins(list(myResults$sims.list$w[ ,t]),
                   x = t,
                   at = t,
                   violin.width = 0.3,
                   col = myCol,
                   add = T,
                   alpha = 0.3,
                   border.col = myCol)
   }#t
   
   params <- dimnames(nimOutput[[1]])[[2]][grep("w",dimnames(nimOutput[[1]])[[2]])]
   for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
   }
   
   ### ====    2.5.phi ====  
   par(mfrow = c(1,1))
   plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "phi")
   axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
   myCol <- c("firebrick3")
   for(t in 1:(nYears-1)){
      plot.violins(list(myResults$sims.list$phi[ ,t]),
                   x = t,
                   at = t,
                   violin.width = 0.2,
                   col = myCol,
                   add = T,
                   alpha = 0.3,
                   border.col = myCol)
   }#t
   
   params <- dimnames(nimOutput[[1]])[[2]][grep("phi",dimnames(nimOutput[[1]])[[2]])]
   for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
   }
   
   ### ====    2.6.p0 ====  
   ## by country and trap-response
   par(mfrow = c(1,2))
   myDev <- c(-0.3, 0.3)
   myCol <- c("blue4", "yellow3")
   
   COUNTIES_AGGREGATEDSubsetsimp <- gSimplify(COUNTIES_AGGREGATEDSubset,tol = 500,topologyPreserve = T)
   COUNTIES_AGGREGATEDSubsetsimp$idunique <- COUNTIES_AGGREGATEDSubset$idunique
   
   for(c in 1:dim(myResults$sims.list$p0)[2]){
      plot(myStudyArea)
      plot(COUNTIES_AGGREGATEDSubset[COUNTIES_AGGREGATEDSubsetsimp$idunique %in% c, ], add=T, col="red")
      
      
      plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.06), type ="n", xaxt="n", xlab = "Years", ylab = "p0")
      axis(1, at = 1:(nYears), labels = years[1:(nYears)])
      
      for(t in 1:nYears){
         plot.violins(list(myResults$sims.list$p0[ ,c,t]),
                      x = t ,
                      at = t ,
                      violin.width = 0.2,
                      col = "red",
                      add = T,
                      alpha = 0.3,
                      border.col = "red")
         # }#g
      }#t
   }#c
   # legend(x = 0, y = 0.25, legend = c("Norway", "Sweden"), fill= myCol)
   
   params <- dimnames(nimOutput[[1]])[[2]][grep("p0", dimnames(nimOutput[[1]])[[2]])[-1]]
   for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
   }
   
   ### ====    2.7.the rest ====
   #beta1 = betaTracks
   #beta2 = betaRoads
   #beta3 = betaSnow
   params <- c("sigma", "dispSigma", "betaResponse", "betaCovs[1]", "betaCovs[2]", "betaCovs[3]", "betaDens")
   for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
   }#i
   dev.off()
}#doall[CM]


## TRY DENSITY ESTIMATES 
gc()
# ORIGIN DENSITY
habbR <- myHabitat$habitat.r#raster::aggregate(myHabitat.list$habitat.r,fact=2)
habbR <- raster::disaggregate(myHabitat$habitat.r, fact=1)
habbR[habbR>0] <- 1
habbdens <- habbR
plot(habbdens)
dim(myResults$sims.list$sxy)
gc()

# RESCALE THE COORDINATES 
sxyScaled <- GridToUTM(data.sxy = myResults$sims.list$sxy,
                       grid.sp = myHabitat$habitat.sp )
gc()

# CHE
plot(habbR)
dim(sxyScaled$data.scaled.xy)
points(sxyScaled$data.scaled.xy[1500,,2,1]~sxyScaled$data.scaled.xy[1500,,1,1], pch=16)
ite <- (dim(sxyScaled$data.scaled.xy)[1]-10000):dim(sxyScaled$data.scaled.xy)[1]
dim(myResults$sims.list$z)


Densi <- list()
for(t in 1:dim(sxyScaled$data.scaled.xy)[4]){
   Densi[[t]] <- EstimateN_v3(habRaster = habbR, # RASTER  FOR WHICH REGION ESTIMATES SHOULD BE ESTIMATED. 
                              posterior.sxy = sxyScaled$data.scaled.xy[,,,t],#sxyScaled$data.scaled.xy[ite,,,t],
                              posterior.z =  myResults$sims.list$z[,,t],#myResults$sims.list$z[ite,,t],
                              alive.states =  c(2),
                              return.all = FALSE,
                              regionEstimates = FALSE)
   print(t)
}
Densi[[t]]$PosteriorsCellsMean

max <- max(unlist(lapply(Densi, function(x) max(x$PosteriorsCellsMean))))
cuts <- seq(0,max,length.out = 100)   #set breaks
col <- rev(terrain.colors(100))

pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_DENSITY.pdf",sep="")))

for(t in 1: dim(sxyScaled$data.scaled.xy)[4]){
   habbdens[] <- Densi[[t]]$PosteriorsCellsMean
   
   plot(habbdens, breaks=cuts, col = col,legend=FALSE, main=years[t]) #p
   plot(myStudyArea.poly,add=T, border=grey(0.5))
   points(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t],],
          pch=16, cex=0.4, col=adjustcolor("black",alpha.f = 0.2))
   plot(habbdens, legend.only=TRUE,breaks=cuts, col=col,
        legend.width = 2,
        axis.args=list(at=round(seq(0, max, length.out = 5),digits = 1),
                       labels=round(seq(0, max, length.out = 5),digits = 1), 
                       cex.axis=0.6),
        legend.args=list(text='Density', side=4, font=2, line=2.5, cex=0.8))
   
}
dev.off()


### ==== 3. PLOT DENSITY MAPS ====
### ====    3.1. PROCESS SXY OUTPUT ====
## List the directories containing bite outputs
outDirectories <- list.files(file.path(myVars$WD, myVars$modelName))[grep("NimbleOut", list.files(file.path(myVars$WD, myVars$modelName)))]
path.list <- file.path(myVars$WD, myVars$modelName, outDirectories)

## Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x)length(list.files(x))))
minBites <- min(numBites)

nimOutput <- list()
for(p in 1:length(path.list)){
   print(path.list[p])
   outfiles <- list.files(path.list[p])
   out <- list()
   for(x in 61:minBites){
      print(x)
      load(file.path(path.list[p], paste("bite_", x, ".RData", sep = "")))
      params.simple <- sapply(strsplit(colnames(mcmcSamples), "\\["), "[", 1)
      parmIndex <- which(params.simple %in% c("sxy","z"))
      out[[x]] <- mcmcSamples[ ,parmIndex] 
   }#x
   out.mx <- do.call(rbind, out)
   nimOutput[[p]] <- as.mcmc(out.mx)
}#p

nimOutput <- as.mcmc.list(nimOutput)
myResults <- ProcessCodaOutput(x = nimOutput)

## RESCALE POSTERIOR sxy
postSxy <- GridToUTM( data.sxy = myResults$sims.list$sxy,
                      grid.sp = myHabitat$habitat.sp)$data.scaled.xy

### ====    3.2. PLOT DENSITY MAPS ====
## NEW GetDensity FUNCTION
GetPopDensityFromPolygon <- function(habPolygon = NULL,
                                     numDiv.x = 1,
                                     numDiv.y = 1,
                                     habResolution = NULL,
                                     habRaster = NULL,
                                     posterior.sxy,
                                     posterior.z,
                                     alive.states,
                                     plot.check = FALSE){
   
   ## Filter out posterior sxy for dead individuals
   deadId <- !which(posterior.z %in% alive.states, arr.ind = TRUE)
   post.x <- posterior.sxy[ , ,1]
   post.x[deadId] <- -999
   posterior.sxy[ , ,1] <- post.x
   
   ## Divide the original polygon using SF package
   habPolygons <- st_make_grid(x = habPolygon, n = c(numDiv.x, numDiv.y))
   habPolygons <- as(habPolygons, "Spatial")
   habPolygons <- habPolygons[!is.na(habPolygons %over% habPolygon)]
   habPolygons <- intersect(habPolygons, as(habPolygon, "SpatialPolygons"))
   
   ## Initialize list of outputs
   thisMean <- thisLowerCI  <- thisUpperCI  <- thisMedian <- list()
   
   for(n in 1:length(habPolygons)){
      
      print(paste("Processing raster bite #", n, " of ", length(habPolygons), sep = ""))
      
      ## Create a raster of the zone of interest
      habRaster <- raster( habPolygons[n],
                           resolution = habResolution,
                           origin = c(extent(habPolygon)[c(1,3)]))
      # res(habRaster) <- habResolution
      habRaster <- rasterize(habPolygons[n], habRaster)
      # origin(habRaster) <- c(150000, 6400000) # c(extent(habPolygon)[c(1,3)])
      
      # Initialize the list of output
      thisMean[[n]] <- thisLowerCI[[n]] <- thisUpperCI[[n]] <- thisMedian[[n]] <- habRaster
      
      if(all(is.na(habRaster[]))){next}
      
      ## Convert habRaster to SpatialGrid object 
      spRaster <- as(habRaster, "SpatialGridDataFrame")
      
      ## Calculate cell-specific density for each iteration
      Density <- apply(X = posterior.sxy, FUN = function(curCoords, inRaster){
         ## Get the cell index for each coordinate
         curCoords <- curCoords[curCoords[ ,1] != -999, ]
         curGridIndeces <- getGridIndex( curCoords, getGridTopology(inRaster),
                                         all.inside = FALSE)
         ## Get rid of coordinate indeces that do not fall in a cell
         curGridIndeces <- curGridIndeces[!is.na(curGridIndeces)]
         ## Create a frequency table of grid indeces
         freqTable <- table(as.character(curGridIndeces))
         ## Initialise an output vector
         outVector <- rep(0, nrow(coordinates(inRaster)))
         outVector[as.integer(names(freqTable))] <- freqTable
         outVector[is.na(spRaster@data$layer)] <- NA
         outVector
      }, inRaster = spRaster, MARGIN = 1)
      
      ## EXTRACT SUMMARY STATISTICS PER HABITAT CELL
      meanDensity <- apply(Density, 1, mean)
      CI <- apply(Density, 1, function(x)quantile(x, c(0.025,0.5,0.975), na.rm = TRUE))
      
      ## FEED IN RASTER
      thisMean[[n]][] <- meanDensity
      thisLowerCI[[n]][] <- CI[1, ]
      thisMedian[[n]][] <- CI[2, ]
      thisUpperCI[[n]][] <- CI[3, ]
   }#n
   
   if(plot.check){
      par(mfrow = c(1,2))
      ## Plot mean density 
      max.cut <- max(meanDensity.r[], na.rm = TRUE)
      cuts <- seq(0, max.cut, length.out = 101) 
      pal <- rev(terrain.colors(length(cuts)))
      plot(meanDensity.r, breaks = cuts, col = pal, main ="Mean density", legend=FALSE, axes=FALSE)
      plot(meanDensity.r, legend.only=TRUE, col=pal,
           legend.width=1, legend.shrink=1,
           axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
                          labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
                          cex.axis=0.6),
           legend.args = list(text=paste("Density (ind.",(habResolution/1000)^2,"km-2)",sep=""), side=4, font=2, line=2.5, cex=0.8))
      plot(habPolygon, add = TRUE)
      
      ## Plot CI width
      max.cut <- max(upCIDensity.r[]-lowCIDensity.r[], na.rm = TRUE)
      cuts <- seq(0, max.cut, length.out = 101) 
      pal <- rev(heat.colors(length(cuts)))
      plot(upCIDensity.r-lowCIDensity.r, breaks = cuts, col = pal, main ="95% CI width", legend = FALSE, axes=FALSE)
      plot(upCIDensity.r, legend.only=TRUE, col=pal,
           legend.width=1, legend.shrink=1,
           axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
                          labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
                          cex.axis=0.6))
      plot(habPolygon, add = TRUE)
   }
   return(list("mean.Density" = test <- do.call(merge, thisMean),
               "median.Density" = do.call(merge,thisMedian),
               "upperCI.Density" = do.call(merge,thisUpperCI),
               "lowerCI.Density" = do.call(merge,thisLowerCI)))
}
## CALCULATE DENSITY
# DENSITY <- list()
# for(t in 1:dim(postSxy)[4]){
#    DENSITY[[t]] <- GetPopDensityFromPolygon(habPolygon = myStudyArea,
#                                             numDiv.x = 10,
#                                             numDiv.y = 20,
#                                             habResolution = 20000,
#                                             posterior.sxy = postSxy[ , , ,t],
#                                             posterior.z = myResults$sims.list$z[,,t],
#                                             alive.states = 2,
#                                             plot.check = F)
#    plot(DENSITY[[t]]$mean.Density)
#    plot(myStudyArea, add = TRUE)
# }#t
# 
# save(DENSITY, file = file.path( myVars$WD,
#                                 myVars$modelName,
#                                 paste("9JM_DENSITY_PD.RData", sep = "")))
## PLOT DENSITY
pdf(file = file.path( myVars$WD,
                      myVars$modelName,
                      paste("9JM_DENSITY_PD.pdf", sep = "")))

max.cut <- max(unlist(lapply(DENSITY, function(x) x$mean.Density[])), na.rm = TRUE)
cuts <- seq(0, max.cut, length.out = 11) 
pal <- rev(terrain.colors(length(cuts)))
# max.cut1 <- max(unlist(lapply(myPopDensity, function(x) x$upperCI.Density[]-x$lowerCI.Density[])), na.rm = TRUE)
# cuts1 <- seq(0, max.cut1, length.out = 101) 
# pal1 <- rev(heat.colors(length(cuts1)))

for(t in 1:nYears){
   par(mfrow = c(1,1))
   
   ## Plot mean density 
   plot(DENSITY[[t]]$mean.Density, breaks = cuts, col = pal,
        main = paste( "N =", round(sum(DENSITY[[t]]$mean.Density[], na.rm = T))),
        legend = F, axes = F)  
   plot(myStudyArea, add = T, border = grey(0.5))
   plot(myHabitat$buffered.habitat.poly, border = grey(0.5), add = TRUE)
   plot(DENSITY[[t]]$mean.Density, legend.only = T, col = pal,
        legend.width = 1, legend.shrink = 1,
        axis.args = list( at = round(seq(0, max.cut, length.out = 11), digits = 2),
                          labels = round(seq(0, max.cut, length.out = 11), digits = 2), 
                          cex.axis = 0.6),
        legend.args = list( text = paste("Density (ind.",(myVars$OUTPUT$mapResolution/1000)^2,"km-2)",sep=""),
                            side = 4, font = 2, line = 2.5, cex = 0.8))
   # ## Plot CI width
   # plot(myPopDensity[[t]]$upperCI.Density-myPopDensity[[t]]$lowerCI.Density, breaks = cuts1, col = pal1, main ="95% CI width", legend = FALSE, axes=FALSE)
   # plot(myStudyArea, add = T, border = grey(0.5))
   # plot(myHabitat$buffered.habitat.poly, border = grey(0.5), add = TRUE)
   # plot(myPopDensity[[t]]$upperCI.Density-myPopDensity[[t]]$lowerCI.Density, legend.only=TRUE, col=pal1,
   #      legend.width=1, legend.shrink=1,
   #      axis.args=list(at=round(seq(0, max.cut1, length.out = 11),digits = 2),
   #                     labels=round(seq(0, max.cut1, length.out = 11),digits = 2), 
   #                     cex.axis=0.6))
   # plot(myHabitat$buffered.habitat.poly, add = TRUE)
}#t
dev.off()

##------------------------------------------------------------------------------