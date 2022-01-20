#' @title Check individual distances between all dectection of an individual and the centroid of its detection
#'
#' @description
#' \code{CheckDistanceDetections} returns a \code{list} object with the distances between detections and centroid of detections and which detection from which id is outside the max.distance if specified 
#'
#' @param y A \code{matrix} or \code{array} of detections 
#'	@param detector.xy A \code{matrix} object with xy coordinates of detectors. 
#' @param max.distance A \code{numeric} distance (radius) used to identify  which detections are too far
#' @param plot.check A \code{logical} plot distance histogram if TRUE
#' @return  a \code{list} object with the distances between detections and centroid of detections and which detection from which id is outside the max.distance if specified 
#' 

CheckDistanceDetections <- function(  y = y
                                    , detector.xy = detector.xy
                                    , max.distance = NULL
                                    , plot.check = TRUE){
   
   ## SET THE DECTECTOR COORDINATES
   detector.xy.sp <- data.frame(detector.xy)
   coordinates(detector.xy.sp) <- detector.xy.sp
   
   ## ADD A TEMPORAL COMPONENT IN Y IF NOT THERE
   if(length(dim(y))!=3){
      y.arr <- array(NA,c(dim(y),1))
      y.arr[,,1] <- y}else{y.arr <- y}
   
   # FIND THE MAX NUMBER OF DETECTIONS 
   max.detections <- apply(y.arr, c(1,3), function(x) sum(x>0))
   
   # FIND WHICH ID WERE DETECTED 
   detected <- apply(max.detections, c(1,2),function(x) as.numeric(x>0))
   max.detections <- max(max.detections)

   # CREATE EMPTY ARRAY TO STORE:
   # DISTANCE BETWEEN DETECTIONS AND ACS
   distance <- array(NA, c(dim(y)[1], max.detections, dim(y.arr)[3]))
   # WHICH DETECTORS HAD A DETECTION 
   detections <- array(NA, c(dim(y)[1], max.detections, dim(y.arr)[3]))
   

   # FOR EACH YEAR 
   for(t in 1:dim(y.arr)[3]){
      detected <- which(detected[,t]>0)
      # FOR EACH ID DETECTED
      for(i in 1:length(detected)){
         # FIND THE DETECTORS WHERE ID WERE DETECTED 
         which.detect <- which(y.arr[detected[i],,t]>0)
         tmp.sp <- detector.xy.sp[which.detect,]
      
         # CONVERT TO SP
         tmp.sxy <- coordinates(tmp.sp)
         tmp.meansxy <- colMeans(tmp.sxy)
         tmp.meansxy <- data.frame(x=tmp.meansxy[1], y=tmp.meansxy[2])
         coordinates(tmp.meansxy) <- tmp.meansxy

         # CALCULATE DISTANCES 
        dist.temp <- gDistance(tmp.sp, tmp.meansxy, byid = T)
        
        # STORE INFORMATION
        # DISTANCE
        distance[detected[i],1:length(dist.temp),t] <-  dist.temp
        # DETECTIONS 
        detections[detected[i],1:length(dist.temp),t] <- which.detect
        }#i
      }#t
   
   # IF PLOT 
   if(plot.check){hist(distance, main="", xlab="Distance from individual centroid detections")}
   
   # PRINT MAXIMUM DISTANCE 
   print(paste("Max distance is", round(max(distance,na.rm=T),digits=2)))
   
   # IF MAX DISTANCE RETURN WHICH ONES ARE OUTSIDE THE MAXIMUM DISTANCE 
   if(!is.null(max.distance)){
      detections.out1 <- detections.out <- which(distance >= max.distance, arr.ind = T)
      print(paste("There are", dim(detections.out1)[1],"detections from",length(unique(detections.out1[,1])), "individuals further away than",max.distance,"units"))
      
      # IF DETECTIONS ARE OUT, TELL WHICH ONES 
      if(dim(detections.out1)[1]!=0){
         for(i in 1:dim(detections.out)[1]){
            detections.out1[i,2] <- detections[detections.out[i,1],detections.out[i,2],detections.out[i,3]]
            }#i
         colnames(detections.out1) <- c("Individual", "detector", "t")
         ID.out <- unique(detections.out1[,1])
         }else{detections.out1 <- ID.out <- NULL}#else
      }else{detections.out1 <- ID.out <- NULL}#else
   
   return(list( detections.out = detections.out1
               , ID.out = ID.out
               , distance = distance))
   }