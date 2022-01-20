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

CheckDistanceDetectionsV3 <- function(  y = y
                                    , detector.xy = detector.xy
                                    , max.distance = NULL
                                    , method = "pairwise"
                                    
                                    , plot.check = TRUE){
   
   
   y<-detector.xy
   
   ##  MAKE CONTAINER FOR FLAGS
   #if(!is.null(Y))
   y.flagged<-y.distance<-list()

  i<-93
      for(i in 1:length(y)){
   
               y.flagged[[i]]<- y.distance[[i]]<-NA #---can't assign NULL in this construct
   
      
   ## SET THE DETECTOR COORDINATES
   detector.sp <- data.frame(detector.xy[[i]])
   dim(detector.sp)
   if(dim(detector.sp)[1]>1){ 
   coordinates(detector.sp) <- detector.sp
detector.index<-1:length(detector.sp)
         this.det<-detector.sp
      dist.det <- gDistance(this.det, this.det, byid = T)
      
      n.too.far<-apply(dist.det,2,function(x){
         sum(x>= max.distance)
      })
      
      to.remove<-NA 
      
      if(any(n.too.far>0)){
         
         
         if(   dim(dist.det)[1]==2 ){
            #print("########")
         to.remove<-2 #---remove detections at the second detector if only 2 detectors and they are too far from each other
      }
      
      if(dim(dist.det)[1]>2 ){
         
   combos<-lapply(1:(length(this.det)-1),function(x){
            combn(1:length(this.det),x)

         })
         
         x<-1
         temp3<-lapply(1:length(combos),function(x) {
            xx<-combos[[x]][,1]
            out<-apply(combos[[x]],2,function(xx){
               this.det2<-this.det[-xx,]
               
               if(!any(gDistance(this.det2, this.det2, byid = T)>=max.distance)){
                  return(xx)
                  #break()
               }
            })
            try(out<-out[,1],silent=TRUE)
            return(out)
         })
         
    
         
         
         to.remove<-temp3[!unlist(lapply(temp3,is.null))][[1]]
         if(is.list(to.remove)){
            to.remove<-to.remove[!unlist(lapply(to.remove,is.null))][[1]]
            #to.remove<-to.remove[[1]]
         }
      }
 
      }
      if(plot.check){
         if(!is.null(to.remove)){
            ext<-data.frame(t(apply(coordinates(this.det),2,min)))
            ext[2,2]<-ext[1,2]+max.distance
            ext[2,1]<-ext[1,1]
            coordinates(ext)<-ext
            plot(this.det,main=i,pch=19,col="grey",cex=1)
            lines(coordinates(ext),main=i,col="green")
            plot(this.det[to.remove,],add=TRUE,pch=19,col="red",cex=0.8)
            box()
         }
      }
      
      
      y.flagged[[i]]<-to.remove

      y.distance[[i]]<-unlist(lapply(to.remove,function(x)max(dist.det[x,])))
      }
      }#--- for-loop (i)


  
   print(paste("Detections removed: ", sum(!is.na(unlist(y.flagged))), " of ",dim(do.call(rbind,y))[1],sep=""))
   print(paste("Individuals affected: ", sum(!unlist(lapply(y.flagged,is.na))), " of ",length(y),sep=""))

      out<-list(y.flagged=y.flagged,y.distance=y.distance) 
      return(out)

   }
