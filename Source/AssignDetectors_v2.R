#' @title Detectors Assignment Function
#'
#' @description
#' \code{AssignDetectorsPD} returns a \code{SpatialPointsDataFrame} object with DNA samples
#'  locations and associated detectors'IDs.
#'
#' @param myData A \code{SpatialPointsDataFrame} object containing the filtered data
#' @param myDetector A \code{list} of SpatialPointsDataFrame objects containing the 
#' detectors locations
#' @param mysubDetectors A \code{list} of SpatialPointsDataFrame objects containing
#'  the subdetectors locations for the binomial models
#'
#' @return A \code{SpatialPointsDataFrame} object with the x and y locations of the
#'  different samples with associated detector's Id 

AssignDetectors_v2 <- function( myData,                              
                              myDetectors,
                              mysubDetectors = NULL){

   ## RETRIEVE NUMBER OF YEARS 
   if(is.null(myData$Year)){myData$Year <- 1}
   years <- min(unique(myData$Year)):max(unique(myData$Year))
   
   ## SET UP THE DETECTORS LIST
   if(!is.list(myDetectors)){myDetectors <- list(myDetectors)}
   
   ## CHECK THAT THE NUMBER OF YEARS MATCH THE LENGTH OF THE DETECTORS LIST
   if(length(myDetectors) != length(years)){
      myDetectors <- lapply(1:length(years), function(x) myDetectors[[1]])
      }
   
   ## IF PAB PROCESS
   if(!is.null(mysubDetectors)){
      ## SET UP THE SUB-DETECTORS LIST
      if(!is.list(mysubDetectors)){mysubDetectors <- list(mysubDetectors)}
      if(length(mysubDetectors) != length(years)){
         mysubDetectors <- lapply(1:length(years), function(x) mysubDetectors[[1]])
      }
      ## INITIALIZE DETECTORS AND SUB-DETECTORS ID COLUMNS
      myData$sub.detector <- NA 
      myData$Detector <- NA
      for(t in 1:length(years)){
         
         #--[RB] using sf package:
         myDetectors.sf<- st_as_sf(mysubDetectors[[t]])
         myData.sf<- st_as_sf(myData[myData$Year == years[t],])
         d<-st_distance(myDetectors.sf, myData.sf)
         myData$Detector[myData$Year == years[t]]  <- apply(d, 2, which.min)            ## Assignate the closest detector to each sample
         
         
         ## CALCULATE DISTANCES BETWEEN SAMPLES AND DETECTORS PER YEAR
         #d <- gDistance(mysubDetectors[[t]], myData[myData$Year == years[t], ], byid = T) 
         ## IDENTIFY THE CLOSEST SUB-DETECTOR FOR EACH SAMPLE
         #temp.sub.det <- apply(d, 1, which.min)  
         temp.sub.det <- apply(d, 2, which.min)  
         ## SET THE ID OF THE SUB-DETECTOR IN THE SAMPLE DATAFRAME
         myData$sub.detector[myData$Year == years[t]] <- temp.sub.det 
         ## SET THE CORRESPONDING DETECTOR ID AS WELL
         main.cell.id <- mysubDetectors[[t]]$main.cell.id[temp.sub.det]     
         myData$Detector[myData$Year == years[t]] <- unlist(lapply(1:length(main.cell.id), function(x) which(myDetectors[[t]]$main.cell.id %in% main.cell.id[x])))
      }#t
   }else{
      myData$Detector <- NA
      for(t in 1:length(years)){

         #--[RB] using sf package:
         myDetectors.sf<- st_as_sf(myDetectors[[t]])
         myData.sf<- st_as_sf(myData[myData$Year == years[t],])
         d<-st_distance(myDetectors.sf, myData.sf)
         myData$Detector[myData$Year == years[t]]  <- apply(d, 2, which.min)            ## Assignate the closest detector to each sample
         
         #d <- gDistance(myDetectors[[t]], myData[myData$Year == years[t], ], byid=T)    ## Calculate the sample/detector distance matrix
         #myData$Detector[myData$Year == years[t]]  <- apply(d, 1, which.min)            ## Assignate the closest detector to each sample
         
         
      }
   }
   
   ## the "Factor" trick...needed to suppress unused factor levels in Id
   myData$Id <- factor(as.character(myData$Id), levels = unique(as.character(myData$Id)))
   
   
   ## Output return
   if(!is.null(mysubDetectors)){                                 # return the number of trials when there's subDetectors 
      n.trials <- lapply(mysubDetectors,function(x)table(x$main.cell.id))
      return(list(myData.sp = myData, n.trials = n.trials))
   } else {
      return(myData)
   }
}
