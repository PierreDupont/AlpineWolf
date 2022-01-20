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

CheckDistanceDetectionsV2 <- function(  y = y
                                    , detector.xy = detector.xy
                                    , max.distance = NULL
                                    , method = "pairwise"
                                    
                                    , plot.check = TRUE){
   
   ##  MAKE CONTAINER FOR FLAGS
   
   y.flagged<-y.distance<-y
   y.flagged[]<-0
   y.distance[]<-NA
   
   ## SET THE DECTECTOR COORDINATES
   detector.sp <- data.frame(detector.xy)
   coordinates(detector.sp) <- detector.sp
   #y
   #dim(y)
   #length(detector.sp)
   detector.index<-apply(y,1,function(x)detector.sp[which(x>0),])
   detections.sp.ls<-     apply(y,1,function(x)detector.sp[which(x>0),])
   
   
      #this.det<-detections.sp.ls[[300]]

   if(method=="pairwise"){
      
      for(i in 1:dim(y)[1]){
         detector.index<-which(y[i,]>0)
         this.det<-detector.sp[detector.index,]
      dist.det <- gDistance(this.det, this.det, byid = T)
      
      n.too.far<-apply(dist.det,2,function(x){
         sum(x>= max.distance)
      })
      
      to.remove<-NULL  
      
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
      
      
      y.flagged[i,detector.index[to.remove]]<-1

      y.distance[i,detector.index[to.remove]]<-unlist(lapply(to.remove,function(x)max(dist.det[x,])))
      }
   }
   
   if(method=="centroid") {
      
      for(i in 1:dim(y)[1]){
         detector.index<-which(y[i,]>0)
         this.det<-detector.sp[detector.index,]
         
         if(length(this.det)>0){
         #TO WEIGH BY NUMBER OF DETECTIONS: centroid<-data.frame(t(data.frame(apply(coordinates(detector.sp[rep(detector.index,y[i,detector.index]),]),2,mean))))
         centroid<-data.frame(t(data.frame(apply(coordinates(detector.sp[detector.index,]),2,mean))))
         
         names(centroid)<-c("x","y")
         coordinates(centroid)<-centroid
         dist.det <- gDistance(this.det, centroid, byid = T)
         dimnames( dist.det )[[1]]<-1:dim(dist.det)[1]
         
         to.remove<-as.numeric(unique(dimnames(dist.det)[[2]][dist.det>max.distance]))
         y.flagged[i,to.remove]<-1
         y.distance[i,to.remove]<-dist.det[,as.character(to.remove)]
         sum( y.flagged[i,]);sum( y[i,]>0)
         }
         }
      
      
      
   }
   
   print(paste("Detections removed: ", sum(y * y.flagged), " of ",sum(y),sep=""))
   print(paste("Individuals affected: ", sum(apply(y.flagged,1,function(x)sum(x)>0)), " of ",sum(apply(y,1,function(x)sum(x)>0)),sep=""))

      out<-list(y.flagged=y.flagged,y.distance=y.distance) 
      return(out)

   }
