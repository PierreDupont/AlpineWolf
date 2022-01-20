#' @title Data set clean-up.
#'
#' @description
#' \code{CleanData} returns a \code{SpatialPointsDataFrame} object with the locations of the DNA samples and associated individual infos.
#'
#' @param dna_samples A \code{DataFrame} object containig raw DNA samples data
#'	@param dead_recoveries A \code{DataFrame} object containing the raw dead recoveries data.
#' @param species_id A \code{Numeric} object containing the id of the focal species (1=Bears;2=Wolverines;3:Wolves).
#' @param country_polygon A \code{SpatialPointsDataFrame} with Country polygon for correct assignment of samples to countries
#' @param threshold_month A \code{Numeric} with initial month of the biological year: 1:January...12=December. all samples with months<threshold get year-1. so they get similar year.  
#' @param keep_dead A \code{logical}  Whether dead recovery should be included (TRUE) or not(FALSE)
#' 
#' @return A \code{SpatialPointsDataFrame} object with the x and y locations of the different samples with associated attributes 
#' Id: the individual id
#' DNAID: the sample Id
#' Year: the sample biological year
#' Date: the actual dinfing date  
#' Species: the identified species
#' Country: the country the sample was found in
#' Birth. alleged birth year (when available; based on dead recoveries)
#' Death: biological year of death (for dead recoveries)
#' Age: alleged age (when available; based on dead recoveries)

CleanData <- function( dna_samples                     ## DNA samples file  
                     , dead_recoveries                 ## Dead recoveries file
                     , species_id = NULL               ## Species code: 1=Bears;2=Wolverines;3=Wolves
                     , country_polygon = NULL          ## Country polygon for correct assignment
                     , threshold_month = 1             ## Initial month of the biological year: 1:January...12=December
                     , keep_dead = TRUE                ## wheather dead recovery should be included or not 
                     , reset_month = FALSE)
   {
   ## Merge DNA and dead recoveries data
   dead_recoveries <- dead_recoveries[ ,c(1,2,3,6,13,14,34,35,44)]
   dna_samples <- dna_samples[ ,c(1,2,7,9,19,20,35,37,38)]
   names(dead_recoveries) <- c("RovBaseId","DNAID","Species","Date","Sex","Age","North","East","Id")
   names(dna_samples) <- c("DNAID","RovBaseId","Origin","Date","North","East","Species","Id","Sex")
   myData <- merge(dead_recoveries, dna_samples, by=c("Id","RovBaseId","DNAID","Species","Sex","Date","East","North"), all=TRUE)
   
   ## Clean the data
   myData <- myData[myData$Id!="",]                                             ## Delete unknown individuals
   myData <- myData[!is.na(myData$Id), ]                                        ## Delete NA individuals
   myData <- myData[!is.na(myData$East), ]                                      ## Delete NA locations
   myData <- myData[!is.na(myData$Date), ]                                      ## Delete NA dates
   myData$Species <- as.numeric(myData$Species)                                 
   if(!is.null(species_id)){myData <- myData[myData$Species==species_id, ]}     ## Delete unwanted species
   
   ## Convert dates to biological years
   myData$Date <- as.POSIXct(strptime(myData$Date, "%d.%m.%Y"))
   myData$Year <- as.numeric(format(myData$Date,"%Y"))
   myData$Month <- as.numeric(format(myData$Date,"%m"))
   myData <- myData[!is.na(myData$Year), ]                                      ## Delete NA dates
   if(!is.null(threshold_month)){
      myData$Year[myData$Month<threshold_month] <- myData$Year[myData$Month<threshold_month] - 1
      if(reset_month){myData$Month[myData$Month<threshold_month] <- myData$Month[myData$Month<threshold_month] + 12}
      }#if 
   
   ## Determine Death and Birth Years
   myData$Age <- as.numeric(as.character(myData$Age))
   myData$RovBaseId <- as.character(myData$RovBaseId)
   myData$Death <- NA
   myData$Death[substr(myData$RovBaseId,1,1)=="M"] <- myData$Year[substr(myData$RovBaseId,1,1)=="M"]
   myData$Birth <- myData$Death-myData$Age
   
   ## Convert samples coordinates to the correct spatial projection
   myData <- SpatialPointsDataFrame(cbind(myData$East, myData$North), myData[,c("Id","Sex","RovBaseId","DNAID","Origin","Species","Date","East","North","Age","Year","Month","Birth","Death")]) 
   proj4string(myData) <- CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
   
   # Overlay with SpatialPolygons to determine the countries 
   if(!is.null(country_polygon)){
      myData$Country <- NA
      myData$Country[!is.na(over(myData, as(country_polygon[which(country_polygon$ISO %in% c("FIN")),],"SpatialPolygons") )) ] <- "F"
      myData$Country[!is.na(over(myData, as(country_polygon[which(country_polygon$ISO %in% c("RUS")),],"SpatialPolygons") )) ] <- "R"
      myData$Country[!is.na(over(myData, as(country_polygon[which(country_polygon$ISO %in% c("GOT")),],"SpatialPolygons") )) ] <- "G"
      myData$Country[!is.na(over(myData, as(country_polygon[which(country_polygon$ISO %in% c("NOR")),],"SpatialPolygons") )) ] <- "N"
      myData$Country[!is.na(over(myData, as(country_polygon[which(country_polygon$ISO %in% c("SWE")),],"SpatialPolygons") )) ] <- "S"
      }#if
   
   # Remove dead recoveries if needed
   if(!keep_dead){myData <- myData[is.na(myData$Death), ]} 
   
   ## the "Factor" trick...needed to suppress unused factor levels in Id
   myData$Id <- factor(as.character(myData$Id), levels = unique(as.character(myData$Id)))
   return(myData)
   }