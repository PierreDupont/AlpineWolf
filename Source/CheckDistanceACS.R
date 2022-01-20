#' @title Check individual distances between consecutive individuals ACS from y detection history
#'
#' @description
#' \code{CheckDistanceDetections} returns a \code{list} object with the distances between detections and centroid of detections and which detection from which id is outside the max.distance if specified 
#'
#' @param y A  \code{array} of detections 
#'	@param detector.xy A \code{matrix} object with xy coordinates of detectors. If detectors changes at every occasions, store each matrix with xy coordinates in a list 
#' @param myStudyArea.sp A \code{SpatialPolygon} of the study area
#' @param plot.check A \code{logical} plot distance histogram/ and movement map if TRUE
#' @return  a \code{list} a matrix of distances between consecutive ACS which id is outside the max.distance if specified 
#' 

CheckDistanceACS <- function(  y = y
                             , detector.sp = detector.sp
                             , myStudyArea.sp = myStudyArea.sp
                             , plot.check = TRUE){
   
   ## N YEARS ##
   n.years <- dim(y)[3]
   ## SET THE DECTECTOR COORDINATES
   
   if(!is.list(detector.sp)){
      detector.sp <- lapply(1:n.years,  function(x) data.frame(detector.sp))
   }
   if(!is.data.frame(detector.sp[[1]])){
      detector.sp <- lapply(1:n.years,  function(x) data.frame(detector.sp[[x]]))
   }

   # FIND WHICH ID WERE DETECTED 
   detected <- apply(y, c(1,3),function(x) as.numeric(sum(x>0)>0))
   max.id.detected <- max(apply(detected,2,function(x) max(which(x==1))),na.rm=T)
   
   # CALCULATE AC CENTER
   AC <- array(NA, c(max.id.detected, 2, dim(y)[3]))
   for( i in 1:max.id.detected){
      for(t in 1:dim(y)[3]){
        if(detected[i,t]==1){
           AC[i,,t] <-  colMeans(detector.sp[[t]][which(y[i,,t]>0),c("x","y")])
        }
      }
   }
   
   # CALCULATE DISTANCE CONSECUTIVE YEARS
   distances <- array(NA, c(max.id.detected, 1, dim(y)[3]-1))
   for( i in 1:max.id.detected){
      for(t in 2:dim(y)[3]){
            distances[i,1,t-1] <-  sqrt((AC[i,1,t] - AC[i,1,t-1])^2 + 
                                       (AC[i,2,t] - AC[i,2,t-1])^2)
               
      }
   }
   
   # CALCULATE DISTANCE CONSECUTIVE YEARS
   occasions.detected <- apply(detected, 1, function(x) which(x!=0) )
   distances.not.consecutive <- array(NA, c(max.id.detected, 1, dim(y)[3]))
   for( i in 1:max.id.detected){
      occ <- occasions.detected[[i]]
      if(length(occ)>1){
         for(t in 2:length(occ)){
         distances.not.consecutive[i,1,occ[t]] <-  sqrt((AC[i,1,occ[t]] - AC[i,1,occ[t-1]])^2 + 
                                        (AC[i,2,occ[t]] - AC[i,2,occ[t-1]])^2)
         
         }
   }
   }
   # plot check 
   if(plot.check){
      # PLOT AC DISTANCES 
      hist(distances.not.consecutive, main = "Distances between consecutive ACS",breaks=40,col=grey(0.5))
      hist(distances, main = "Distances between consecutive ACS",breaks=40,
           col=adjustcolor(grey(0.11),alpha.f = 0.8),add=T)
      legend("topright", fill=c(grey(0.5),grey(0.11)),cex=0.7,
             legend = c("Distances  not consecutive ACS","Distancesconsecutive ACS" ))
      # PLOT MOVEMENT 
      par(mar=c(1,1,1,1))
      plot(myStudyArea.sp, border="white", col=grey(0.8))
      col <- rainbow(max.id.detected)
      for(i in 1:max.id.detected){
         for(t in 1:dim(y)[3]){
         points(AC[i,2,t]~AC[i,1,t], bg=col[i], pch=21, cex=0.8,col="black" )
         }
         
         for(t in 2:dim(y)[3]){
            suppressWarnings(arrows(x0 = AC[i,1,t-1], x1 = AC[i,1,t],
                   y0 = AC[i,2,t-1], y1= AC[i,2,t],
                   col=col[i], lwd = 1, length = 0.07, angle = 30))
         }
      }
   }

   return(list(distances.consecutive = distances,
          distances.not.consecutive=distances.not.consecutive))
}