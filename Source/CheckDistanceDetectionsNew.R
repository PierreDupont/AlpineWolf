#' @title Check individual distances between all detections of an individual and the centroid of its detection
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

CheckDistanceDetectionsNew <- function(  y = y
                                    , detector.xy = detector.xy
                                    , max.distance = NULL
                                    , plot.check = TRUE){
   
   ## SET THE DECTECTOR COORDINATES
   detector.xy.sp <- data.frame(detector.xy)
   coordinates(detector.xy.sp) <- detector.xy.sp
   
   ## ADD A TEMPORAL COMPONENT IN Y IF NOT THERE
   if(length(dim(y))!=3){
      y.arr <- array(NA,c(dim(y),1))
      y.arr[,,1] <- y}else{y.arr <- y
      }
   
   # FIND THE MAX NUMBER OF DETECTIONS 
   max.detections <- apply(y.arr, c(1,3), function(x) sum(x>0))
   
   # FIND WHICH ID WERE DETECTED 
   detected.orig <- apply(max.detections, c(1,2),function(x) as.numeric(x>0))
   max.detections <- max(max.detections)

   # CREATE EMPTY ARRAY TO STORE:
   # DISTANCE BETWEEN DETECTIONS AND ACS
   distance <- array(NA, c(dim(y)[1], max.detections, dim(y.arr)[3]))
   # WHICH DETECTORS HAD A DETECTION 
   detections <- array(NA, c(dim(y)[1], max.detections, dim(y.arr)[3]))
   # DETECTIONS FLAGGED FOR REMOVAL
   flagged.removal <- array(FALSE, c(dim(y)[1], max.detections, dim(y.arr)[3]))
   
   if(FALSE){
   #############################################
   ##-----BASED ON DISTANCE FROM THE CENTROID
   #############################################
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
      
}#--TOGGLE
   
   if(TRUE){
   #############################################
   ##-----SEQUENTIAL REMOVAL OF OUTERLYING POINTS
   #############################################
   # FOR EACH YEAR 
   for(t in 1:dim(y.arr)[3]){
      detected <- which(  detected.orig [,t]>0)
      # FOR EACH ID DETECTED
      #par(mfrow=c(3,3),mar=c(1,1,1,1))
      for(i in 1:length(detected)){
         # FIND THE DETECTORS WHERE ID WERE DETECTED 
         which.detect <- which(y.arr[detected[i],,t]>0)
         tmp.sp <- detector.xy.sp[which.detect,]
      

         # CONVERT TO SP
         tmp.sxy <- coordinates(tmp.sp)
         tmp.meansxy <- colMeans(tmp.sxy)
         tmp.meansxy <- data.frame(x=tmp.meansxy[1], y=tmp.meansxy[2])
         coordinates(tmp.meansxy) <- tmp.meansxy

         proj4string(tmp.sp)<-"+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
         
         # CALCULATE DISTANCES 
        #dist.temp <- gDistance(tmp.sp, tmp.meansxy, byid = T)
        dist.temp <- gDistance(tmp.sp, tmp.sp, byid = T)
        
        n.too.far<-apply(dist.temp,2,function(x){
           sum(x>= max.distance)
         })
        
        to.remove<-NULL  
        
        if(any(n.too.far>0)){
        
        #combn(1:length(tmp.sp),length(tmp.sp)-1)#, m, FUN = NULL, simplify = TRUE, ...)
        combos<-lapply(1:(length(tmp.sp)-1),function(x){
           combn(1:length(tmp.sp),x)
                     #tmp.sp[-]
            })
        
        x<-1
        temp3<-lapply(1:length(combos),function(x) {
           xx<-combos[[x]][,1]
           out<-apply(combos[[x]],2,function(xx){
              tmp.sp2<-tmp.sp[-xx,]
              
              # ext<-data.frame(t(apply(coordinates(tmp.sp2),2,min)))
              # ext[2,]<-ext[1,]+max.distance
              # coordinates(ext)<-ext
              # plot(coordinates(ext))
              # plot(tmp.sp2,main=i,pch=19,col="grey",cex=1,add=TRUE)
              # points(coordinates(ext),main=i,col="green",type="b")
              # box()
             
              if(!any(gDistance(tmp.sp2, tmp.sp2, byid = T)>=max.distance)){
                 return(xx)
                 #break()
                 }
                 })
          try(out<-out[,1],silent=TRUE)
           return(out)
                })

to.remove<-unlist(temp3[!unlist(lapply(temp3,is.null))][[1]][[1]])
}
if(plot.check){
        if(!is.null(to.remove)){
   ext<-data.frame(t(apply(coordinates(tmp.sp),2,min)))
   ext[2,2]<-ext[1,2]+max.distance
   ext[2,1]<-ext[1,1]
   coordinates(ext)<-ext
   plot(tmp.sp,main=i,pch=19,col="grey",cex=1)
   lines(coordinates(ext),main=i,col="green")
   plot(tmp.sp[to.remove,],add=TRUE,pch=19,col="red",cex=0.8)
   box()
        }
}
        #lapply(combos,function(x))
        
        # STORE INFORMATION
        # DISTANCE
        distance[detected[i],1:dim(dist.temp)[2],t] <-  apply(dist.temp,2,max)
        # DETECTIONS 
        detections[detected[i],1:dim(dist.temp)[2],t] <- which.detect
        
        # DETECTIONS FLAGGED FOR REMOVAL
        flagged.removal[detected[i],to.remove,t] <- TRUE
        }#i
      }#t
   
      
      # IF PLOT 
      if(plot.check){hist(distance, main="", xlab="Distance from individual centroid detections")}
      
      # PRINT MAXIMUM DISTANCE 
      print(paste("Max distance is", round(max(distance,na.rm=T),digits=2)))
      
      # IF MAX DISTANCE RETURN WHICH ONES ARE OUTSIDE THE MAXIMUM DISTANCE 
      if(!is.null(max.distance)){
         
         
         detections.out1 <- detections.out <- which(flagged.removal, arr.ind = T)
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
      
   }#----TOGGLE    

    
  #points( detector.xy [detections.out1[detections.out1[,"Individual"]==503,"detector"],],add=TRUE)
   
   return(list( detections.out = detections.out1
               , ID.out = ID.out
               , distance = distance))
   }
