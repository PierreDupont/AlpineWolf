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

CleanData_v2 <- function( dna_samples                     ## DNA samples file  
                        , dead_recoveries                 ## Dead recoveries file
                        , species_id = NULL               ## Species code: 1=Bears;2=Wolverines;3=Wolves
                        , country_polygon = NULL          ## Country polygon for correct assignment
                        , threshold_month = 1             ## Initial month of the biological year: 1:January...12=December
                        , keep_dead = TRUE                ## wheather dead recovery should be included or not 
                        , reset_month = FALSE
                        , keep_unidentified_id = FALSE
                        , age.label.lookup = NULL)
   {
   
   ## Merge DNA and dead recoveries data## Rename Dead Recoveries
   names(dead_recoveries)[grep(pattern = "Individ",x = names(dead_recoveries))] <- "Id"
   names(dead_recoveries)[grep(pattern = "RovbaseID",x = names(dead_recoveries),fixed = TRUE)] <- "RovBaseId"
   names(dead_recoveries)[grep(pattern = "DNAID..Proeve",x = names(dead_recoveries),fixed = TRUE)] <- "DNAID"
   names(dead_recoveries)[grep(pattern = "Art",x = names(dead_recoveries))] <- "Species"
   names(dead_recoveries)[grep(pattern = "Kjoenn",x = names(dead_recoveries))] <- "Sex"
   names(dead_recoveries)[grep(pattern = "Doedsdato",x = names(dead_recoveries))] <- "Date"
   names(dead_recoveries)[grep(pattern = "Oest..UTM33",x = names(dead_recoveries))] <- "East"
   names(dead_recoveries)[grep(pattern = "Nord..UTM33",x = names(dead_recoveries))] <- "North"
   names(dead_recoveries)[grep(pattern = "Alder",x = names(dead_recoveries))] <- "Age"
   names(dead_recoveries)[grep(pattern = "Doedsarsak..annet",x = names(dead_recoveries))] <- "DeathCause_2"
   names(dead_recoveries)[grep(pattern = "Doedsarsak",x = names(dead_recoveries))] <- "DeathCause"
   
   ## Rename DNA Samples
   names(dna_samples)[grep(pattern = "Individ",x = names(dna_samples))] <- "Id"
   names(dna_samples)[grep(pattern = "RovbaseID..Proeve",x = names(dna_samples))] <- "RovBaseId"
   names(dna_samples)[grep(pattern = "DNAID..Proeve",x = names(dna_samples))] <- "DNAID"
   names(dna_samples)[grep(pattern = "Art..Analyse",x = names(dna_samples))] <- "Species"
   names(dna_samples)[grep(pattern = "Kjoenn",x = names(dna_samples))] <- "Sex"
   names(dna_samples)[grep(pattern = "Funnetdato",x = names(dna_samples))] <- "Date"
   names(dna_samples)[grep(pattern = "Oest..UTM33",x = names(dna_samples))] <- "East"
   names(dna_samples)[grep(pattern = "Nord..UTM33",x = names(dna_samples))] <- "North"
   names(dna_samples)[grep(pattern = "Alder",x = names(dna_samples))] <- "Age"
   names(dna_samples)[grep(pattern = "DNA.proeveleverandoer",x = names(dna_samples))] <- "Origin"
   
   #myData <- merge(dead_recoveries, dna_samples, by=c("Id","RovBaseId","DNAID","Species","Sex","Date","East","North"), all = TRUE)
   myData <- dna_samples
   myData$Death<-NA
   ## Clean the data
   if(keep_unidentified_id==FALSE){                                             ## MAY want to keep all observation, even the unidentified IDS
      myData <- myData[myData$Id!="",]                                          ## Delete unknown individuals
      myData <- myData[!is.na(myData$Id), ]                                     ## Delete NA individuals
   }
   myData <- myData[!is.na(myData$East), ]                                      ## Delete NA locations
   myData <- myData[!is.na(myData$Date), ]                                      ## Delete NA dates
   if(!is.null(species_id)){myData <- myData[myData$Species %in% species_id, ]} ## Delete unwanted species
   
   ## Convert dates to biological years
   myData$Date <- as.POSIXct(strptime(myData$Date, "%d.%m.%Y"))
   myData$Year <- as.numeric(format(myData$Date,"%Y"))
   myData$Month <- as.numeric(format(myData$Date,"%m"))
   myData <- myData[!is.na(myData$Year), ]                                      ## Delete NA dates
   if(!is.null(threshold_month)){
      myData$Year[myData$Month<threshold_month] <- myData$Year[myData$Month<threshold_month] - 1
      if(reset_month){myData$Month[myData$Month < threshold_month] <- myData$Month[myData$Month < threshold_month] + 12}
      }#if 
   
   ## Determine Death and Birth Years
   # myData$Age <- suppressWarnings(as.numeric(as.character(myData$Age))) 
   # myData$RovBaseId <- as.character(myData$RovBaseId)
   # myData$Death <- NA
   # myData$Death[substr(myData$RovBaseId,1,1)=="M"] <- myData$Year[substr(myData$RovBaseId,1,1)=="M"]
   # myData$Birth <- myData$Death-myData$Age
   # 
   # ## Reconstruct minimal & maximal ages
   # myData$Age.orig <- myData$Age
   # if(!is.null(age.label.lookup)){
   #    temp <- temp1 <- as.character(levels(myData$Age.orig))  ## list age levels
   #    temp <- toupper(temp)                                   ## Upper case all
   #    temp <- gsub("\\s", "", temp)                           ## Remove blank spaces
   #    myData$Age.orig2 <- myData$Age.orig
   #    levels(myData$Age.orig2) <- temp
   #    myData <- merge( myData, age.label.lookup[ ,-1],
   #                     by.x = "Age.orig2",
   #                     by.y = "age.label",
   #                     all.x = TRUE)                          ## Merge with info from lookup table
   #    
   #    ## FILL IN THE REST OF THE AGES FROM FOR NUMERIC RECORDS
   #    numeric.age.records <- which(!is.na(as.numeric(as.character(myData$Age.orig2))) & !is.na(myData$Age.orig2))
   #    myData[numeric.age.records, c("min.age","max.age","age")] <- floor(as.numeric(as.character(myData$Age.orig2[numeric.age.records])))
   # }
   
   ## Convert samples coordinates to the correct spatial projection
   myData <- SpatialPointsDataFrame(coords = cbind(myData$East, myData$North), data = myData) 
   proj4string(myData) <- proj4string(country_polygon)#CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
   
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