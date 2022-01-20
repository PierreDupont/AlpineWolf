#' @title Function to set initial values of XY coordinates of ACS
#'
#' @description
#' \code{MakeInitsXY} returns a matrix object with with the x coordinates ([,1]) and y coordinates  ([,2]). it returns the location of a detection for a detected individual and a random location for augmented individuals.
#' for individuals that remained undetected during one or several years, the detections from before/after are used as the centroid of the buffer (radius determinde by the dist move).
#' 
#' @param y \code{matrix} or \code{array}  with individual detections. row=Ids, col= detectors and if arrray=  [,,t] time. 
#' @param detector.xy \code{matrix} or \code{array} with coordinates of detectors. if array, the 3rd dimension corresponds to years
#' @param habitat.r \code{raster} with the raster of the habitat (0/1)
#' @param dist.move \code{numeric} with a radius value of how much an individual can move from year to year. sxy t+1 will be drawn from a buffer centered on average sxy t
#' @param ydead \code{array}  with individual dead recoveries row=Ids, col= detectors and if arrray=  [,,t] time. 
#' @param detector.xyDead \code{matrix} or \code{array} with coordinates of detectors. if array, the 3rd dimension corresponds to years
#'
#'
#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' MakeInitsXY(   y = y
#'              , detector.xy = data$detector.xy
#'              , habitat.r = habitat.r
#'              , dist.move = dist.move )

MakeInitsXY_v2 <- function(    y = y
                           , 
                           detector.xy = detector.xy
                           ,
                           habitat.r =  myHabitat$habitat.r
                           ,
                           dist.move = dist.move
                           ,
                           ydead = NULL
                           ,
                           detector.xyDead = NULL
                           ,
                           seed = FALSE
                           ,
                           buffer.r = 10000
){
   
   ##PREPARE INPUT  
   n.years <- ifelse(length(dim(y))>=3, dim(y)[3], 1)
   n.individuals <- dim(y)[1]
   n.detectors <- dim(detector.xy)[1]
   if(length(dim(y)) == 2){y <- array(y, c(n.individuals, n.detectors, n.years))}
   if(length(dim(detector.xy)) == 2){detector.xy <- array(detector.xy, c(n.detectors, 2, n.years))}
   sxy <- array(NA, c(n.individuals, 2, n.years))
   
   # CHECK WHO/WHEN IDS ARE DETECTED 
   detected<-NULL
   try(
   detected <- apply(y, c(1,3), function(x) any(x >= 1))
   ,silent=TRUE)
   #---[RB]: dealing with size alocation issues (breaking things down into smaller bite sizes)
   
   if(is.null(detected)){
   temp<-round(seq(1,dim(y)[1],length.out=10))
   temp[length(temp)]<-dim(y)[1]
   temp<-cbind(temp[-length(temp)],temp[-1])
   temp[2:dim(temp)[1],1]<-temp[2:dim(temp)[1],1]+1
   detected<-lapply(1:dim(temp)[1],function(i){
   detected <- apply(y[temp[i,1]:temp[i,2],,], c(1,3), function(x) any(x >= 1))
   })
   detected<-do.call(rbind,detected)
   }
   
   id.detec <- which(apply(detected, 1, function(x) sum(x==FALSE)!=n.years))
   
   if(!is.null(ydead)){
      n.detectors.dead <- dim(detector.xyDead)[1]
      if(length(dim(detector.xyDead)) == 2){detector.xyDead <- array(detector.xyDead, c(n.detectors.dead, 2, n.years))}
      detected.dead <-  ydead > 0
      id.dead <- which(apply(ydead,1, function(x) any(x > 0)) ==T)
      id.detec <- unique(c(id.detec, id.dead ))
   }
   
   # CREATE THE HABITAT POLYGON
   buffered.habitat.poly <- aggregate(rasterToPolygons(habitat.r, fun=function(x){x>0}))
   crs <- CRS(proj4string(buffered.habitat.poly))
   
   ## DETECTED INDIVIDUALS 
   for(i in id.detec){
      
      # IF AT LEAST ONE DETECTION ALIVE
      if(sum(detected[i,])>0){
         # YEAR WITH DETECTIONS 
         detected.years <- which(detected[i,]==TRUE)
         #GET AVERAGE DETECTIONS
         for(t in detected.years){
            tmpsxy <- detector.xy[which(y[i,,t]>0),,t] 
            sxy[i,,t] <- if(is.null(dim(tmpsxy))){
               tmpsxy }else{ colMeans(tmpsxy)}  
            
            ## if AVERAGE DETECTION HAPPENS TO BE OUTSIDE THE HABTIAT, THEN FIND A CLOSEST PLACE AROUND
            if(habitat.r[cellFromXY(habitat.r, sxy[i,,t])] ==0){
               sp <- SpatialPoints(rbind(sxy[i,,t]), proj4string = crs)
               buffer <- gBuffer(sp, width = dist.move+buffer.r, byid = T )
               int <- raster::intersect(buffer, buffered.habitat.poly)
               if(seed){
                  set.seed(i + seed)
               }
               sxy[i,,t] <- t(coordinates(spsample(int, n=1, type="random", iter=500)))
               
            }
            
            
         }
         # YEAR WITH NODETECTIONS
         if(sum(detected[i,])!=n.years){
             notdetected.years <- which(detected[i,]==FALSE)
            
             for(t in notdetected.years){
                if(t < min(detected.years)){
                   sxy[i,,t] <- sxy[i,,min(detected.years)] 
                }
                if(t > max(detected.years)){
                   sxy[i,,t] <- sxy[i,,max(detected.years)]  
                }
                if(t < max(detected.years) & t > min(detected.years) ){
                      sxy[i,,t] <- sxy[i,,t-1]  
                }
                
                 sp <- SpatialPoints(rbind(sxy[i,,t]), proj4string = crs)
                 buffer <- gBuffer(sp, width = dist.move, byid = T )
                 int <- raster::intersect(buffer, buffered.habitat.poly)
                 if(seed){
                    set.seed(i + seed)
                 }
                 sxy[i,,t] <- t(coordinates(spsample(int, n=1, type="random", iter=500)))
                 
             }
             

            #CALCULATE AVERAGE DETECTIONS OF ALL YEARS 
            # notdetected.years <- which(detected[i,]==FALSE)
            # averagesxy <- sxy[i,,detected.years] 
            # averagesxy <-  if(is.null(dim(averagesxy))){
            #    averagesxy }else{ rowMeans(averagesxy)}  
            # DRAW A BUFFER
            # sp <- SpatialPoints(rbind(averagesxy), proj4string = crs)
            # buffer <- gBuffer(sp, width = dist.move )
            # int <- intersect(buffer, buffered.habitat.poly)
            # DRAW SXY
         
         }
      }
      
      ## IF DEAD
      if(!is.null(ydead)){
         # IF AT LEAST ONE DETECTION DEAD
         if(sum(detected.dead[i,])>0){
            # YEAR WITH DETECTIONS 
            detected.years <- which(detected.dead[i,]==TRUE)
            #GET SXY DEAD 
            for(t in detected.years){
               sxy[i,,t] <- detector.xyDead[ydead[i,t],,t] 
            }
           # AFTER BEING DEAD, STAY WHERE IT IS DEAD
            detected.deadi <- min(which(detected.dead[i,]==TRUE))
           if(detected.years < n.years){
              #make a buffer 
              sp <- SpatialPoints(rbind(sxy[i,,detected.deadi]), proj4string = crs)
              buffer <- gBuffer(sp, width = dist.move )
              int <- raster::intersect(buffer, buffered.habitat.poly)
              if(seed){
                 set.seed(i + seed)
              }
              # DRAW SXY
              sxy[i,,detected.deadi:n.years] <- t(coordinates(spsample(int, n=length(detected.deadi:n.years), 
                                                                       type="random", iter=500)))
              }
            
           
            
            #IF INDIVIDUAL ONLY DETECTED DEAD 
            if(sum(detected[i,])==0){
               #CALCULATE AVERAGE DETECTIONS OF ALL YEARS 
               notdetected.years <- which((detected[i,]+detected.dead[i,])==0)
               averagesxy <- sxy[i,,detected.years] 
               # DRAW A BUFFER
               sp <- SpatialPoints(rbind(averagesxy), proj4string = crs)
               buffer <- gBuffer(sp, width = dist.move )
               int <- raster::intersect(buffer, buffered.habitat.poly)
               if(seed){
                  set.seed(i + seed)
               }
               # DRAW SXY
               sxy[i,,notdetected.years] <- t(coordinates(spsample(int, n=length(notdetected.years), type="random", iter=500)))
            }
         }
      }
   }
   
   #AUGMENTED IDS 
   id.notdetec <- c(1:n.individuals)[(c(1:n.individuals) %in% id.detec)==FALSE]
   if(length(id.notdetec)!=0){
   start.loca <- spsample(buffered.habitat.poly, n=length(id.notdetec),type="random",iter=500)
   # DRAW A BUFFER
   buffer.augm <- gBuffer(start.loca,width = dist.move, byid = T)
   buffered.habitat.poly <- gBuffer(buffered.habitat.poly, width = -1)
   int.buff <- raster::intersect(buffered.habitat.poly, buffer.augm)
   
   for(i in 1:length(id.notdetec)){
      #DRAW COORDINATES FROM WITHIN THE BUFFER 
      if(seed){
         set.seed(i + seed)
      }
      coords <- spsample(int.buff[i,], n=n.years, type="random",iter=500)
      sxy[id.notdetec[i],,] <- t(coordinates(coords))
      #model$habitat.mx[trunc(thisSXY[2])+1, trunc(thisSXY[1])+1]    ## coming from habitat array
      
   }
   }
   # ADD SOME RANDOM NOISE TO AVOID NON-ZERO MOVEMENT
   for(t in 1:n.years){
      if(seed){
         set.seed(t + seed)
      }
      sxy[,2,t] <-  sxy[,2,t] +jitter(0,factor = 0.01)
   }
   return(sxy)
   
}



