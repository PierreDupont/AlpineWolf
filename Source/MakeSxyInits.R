#' @title Check individual distances between all dectection of an individual 
#' and the centroid of its detection
#'
#' @description
#' \code{CheckDistanceDetections} returns a \code{list} object with the distances 
#' between detections and centroid of detections and which detection from which id 
#' is outside the max.distance if specified 
#'
#' @param y A \code{matrix} or \code{array} of detections 
#'	@param detector.xy A \code{matrix} object with xy coordinates of detectors. 
#' @param max.distance A \code{numeric} distance (radius) used to identify  which 
#' detections are too far
#' @param plot.check A \code{logical} plot distance histogram if TRUE
#' @return  a \code{list} object with the distances between detections and centroid
#' of detections and which detection from which id is outside the max.distance if 
#' specified 

MakeSxyInits <- function( y, 
                          habitatPoly,
                          detCoords = NULL,
                          detNum = NULL,
                          maxDetDist = NULL,
                          maxMoveDist = NULL,
                          y.dead = NULL,
                          detCoords.dead = NULL,
                          detNum.dead = NULL,
                          plot.check = TRUE){
   
   ## ---- I. PREPARE INPUT ----
   if(!is.null(detCoords)){ 
      if(length(dim(y)) == 2){
         y <- array(y, c(dim(y), 1))
         if(!is.null(y.dead)){
            y.dead <- matrix(y.dead, length(y.dead), 1)
         }
      }
      nYears <- dim(y)[3]
      nInd <- dim(y)[1]
      detCoords <- array( detCoords,
                          c(dim(detCoords)[1], dim(detCoords)[2], nYears))
   }
   
   if(!is.null(detNum)){
      if(length(dim(y)) == 3){
         y <- array(y, c(dim(y), 1))
         if(!is.null(y.dead)){
            y.dead <- array(y.dead, c(dim(y.dead), 1))
         }
      }
      nYears <- dim(y)[4]
      nInd <- dim(y)[3]
      detNum <- matrix(detNum, nInd, nYears)
   }
   
   ## ---- II. PREPARE OUTPUT OBJECTS ----
   sxy <- array(NA, c(nInd, 2, nYears))
   detDistList <- moveDistList <- NULL
   
   ## ---- Loop over individuals ----
   for(i in 1:nInd){
      ## Initialize the list of polygons of yearly possible sxy locations
      potentialSxy <- detSpatialPoints <- list()
      if(plot.check){plot(habitatPoly, col = "gray80")}
      
      ## ---- III. FIRST PASS ----
      ## CREATE POLYGONS OF POTENTIAL SXY LOCATIONS BASED ON YEARLY DETECTIONS 
      for(t in 1:nYears){
         ## Initialize the polygon of possible sxy locations
         potentialSxy[[t]] <- habitatPoly
         detSpatialPoints[[t]] <- NA
         temp <- NULL
         
         ## Retrieve detection coordinates (if any)
         if(!is.null(detCoords)){
            if(any(y[i,,t] > 0)){temp <- cbind.data.frame(detCoords[which(y[i, ,t] > 0),1,t],
                                                          detCoords[which(y[i, ,t] > 0),2,t])}
            if(!is.null(y.dead)){
               detCoords.dead <- array(detCoords.dead,
                                       c(dim(detCoords.dead)[1],
                                         dim(detCoords.dead)[2],
                                         nYears))
               if(y.dead[i,t] > 0){temp <- cbind.data.frame(detCoords.dead[y.dead[i,t],1,t],
                                                            detCoords.dead[y.dead[i,t],2,t])}
               if(any(y.alive[i,,t] > 0) & y.dead[i,t] > 0){
                  stop(paste("Individual",i,"was detected ALIVE & DEAD at time",t))}
            }
         }#if
         
         ## Retrieve detection coordinates (if any)
         if(!is.null(detNum)){
            if(detNum[i,t] > 0){temp <- cbind.data.frame(y[1:detNum[i,t],1,i,t],
                                                         y[1:detNum[i,t],2,i,t])}
            if(!is.null(y.dead)){
               detNum.dead <- matrix(detNum.dead, nInd, nYears)
               if(detNum.dead[i,t] > 0){temp <- cbind.data.frame(y.dead[1,i,t], 
                                                                 y.dead[2,i,t])}
               if(detNum[i,t] > 0 & detNum.dead[i,t] > 0){
                  stop(paste("Individual",i,"was detected ALIVE & DEAD at time",t))}
            }
         }#if
         
         ## If multiple detections  
         if(!is.null(temp) & !is.null(maxDetDist)){
            
            ## Transform detection coordinates to spatial points
            detSpatialPoints[[t]] <- SpatialPoints( temp,
                                                    proj4string = CRS(projection(habitatPoly)))
            
            ## Loop over detections
            for(n in 1:length(detSpatialPoints[[t]])){
               ## Buffer the current detection w/ maximum distance allowed 
               ## between sxy and detections
               bufferedPoint <- gBuffer( spgeom = detSpatialPoints[[t]][n],
                                         width = maxDetDist)
               if(plot.check){
                  plot(bufferedPoint, col = rgb(0.41,0.51,0.55,0.2), add = T)
                  plot(detSpatialPoints[[t]][n], add = TRUE)
               }
               
               ## Determine the intersection between the buffered detection 
               ## and the current polygon of possible sxy locations
               tempSxy <- crop(habitatPoly, extent(bufferedPoint))
               tempSxy <- gIntersection(bufferedPoint, tempSxy)
               ## Check that an intersection exists
               if(is.null(tempSxy)){
                  stop(paste("ERROR: maxDetDist too small to allow for the
                             detections of individual ", i, " at time ", t,
                             "!", sep = ""))
               }#if
               potentialSxy[[t]] <- tempSxy
               ## Check that there is no weird intersection problems
               if(class(potentialSxy[[t]])=="SpatialCollections"){
                  potentialSxy[[t]] <- potentialSxy[[t]]@polyobj
                  }
            }#n
         }#if
      }#t
      
      ## ---- IV. SECOND PASS ---- 
      ## CREATE POLYGONS OF POTENTIAL SXY LOCATIONS BASED ON 
      ## THE DISTANCES BETWEEN OTHER POTENTIAL YEARLY SXY LOCATIONS
      for(t in 1:nYears){
         if(nYears > 1 & !is.null(maxMoveDist)){
            otherYears <- (1:nYears)[-t]
            for(t2 in otherYears){
               if(class(potentialSxy[[t2]]) == "SpatialPoints" | class(detSpatialPoints[[t2]]) == "SpatialPoints"){
                  ## Buffer polygons of potential sxy locations w/ maximum distance 
                  ## allowed between two consecutive sxy locations multiplied by 
                  ## the number of time intervals
                  bufferedPoly <- rgeos::gBuffer( spgeom = potentialSxy[[t2]],
                                                  width = maxMoveDist * abs(t-t2))
                  if(class(detSpatialPoints[[t2]]) == "SpatialPoints"){
                     bufferedPoly2 <- rgeos::gBuffer( spgeom = detSpatialPoints[[t2]],
                                                      width = maxMoveDist * abs(t-t2))
                     bufferedPoly <- rgeos::gIntersection(bufferedPoly, bufferedPoly2)
                  }#if
                  if(plot.check){plot(bufferedPoly, col = rgb(0,1,0,0.1), add = T)}
                  ## Intersect polygons of buffered potential sxy locations with 
                  ## current polygon of potential sxy locations
                  tempSxy <- rgeos::gIntersection(bufferedPoly, potentialSxy[[t]])
                  ## Check that an intersection exists
                  if(is.null(tempSxy)){stop(paste("ERROR: maxMoveDist too small 
                                                  to allow for the detections of 
                                                  individual ", i, " at time ", t,
                                                  "!", sep = ""))}#if
                  potentialSxy[[t]] <- tempSxy
                  ## Check that there is no weird intersection problems
                  if(class(potentialSxy[[t]])=="SpatialCollections"){
                     potentialSxy[[t]] <- potentialSxy[[t]]@polyobj
                     }
               }#if
            }#t2
         }#if
         
         ## Plot the polygon of possible sxy locations
         if(plot.check){plot(potentialSxy[[t]], col = rgb(1,0,0,0.5), add = T)}
         ## Sample one location from the polygon of potential locations
         potentialSxy[[t]] <- spsample( x = potentialSxy[[t]],
                                        n =  1,
                                        type = "random",
                                        iter = 100)
         ## Plot the proposed sxy location
         if(plot.check){plot( potentialSxy[[t]],
                              col = "blue", pch = 19,
                              cex = 1.3, add = T)}
         ## Store the proposed sxy
         sxy[i, ,t] <- coordinates(potentialSxy[[t]])
      }#t
   }#i
   return(sxy)
}


