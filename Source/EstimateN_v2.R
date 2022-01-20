EstimateN_v2 <- function(habPolygon = NULL,
                              habResolution = NULL,
                              habRaster = NULL,
                              posterior.sxy,
                              posterior.z,
                              alive.states,
                              plot.check=FALSE,
                         return.all=FALSE
                         ){
   
   # Retrieve the resolution of the raster
   if(!is.null(habRaster)){
      habResolution <- res(habRaster)
      }
   
   # Create a raster of the zone of interest
   if(is.null(habRaster)){
      habRaster <- raster(extent(habPolygon))
      res(habRaster) <- habResolution
      habRaster <- rasterize(habPolygon, habRaster)
      }
   

      # sxy.<-lapply(1:dim(posterior.sxy)[1],function(i)posterior.sxy[i,,])
      # sxy.<-do.call(rbind,sxy.)
      # dim(sxy.)
      # 
      # z.<-lapply(1:dim(posterior.z)[1],function(i)posterior.z[i,])
      # z.<-do.call(c,z.)
      # dim(z.)
      # 
      # keep<-z.%in%alive.states
      # z.<-z.[keep]
      # sxy.<-sxy.[keep,]
      # 
      # sxy<-data.frame(sxy.[1:100000,])#data.frame(posterior.sxy[i,posterior.z[i,]%in%alive.states,])
      # 
      # coordinates(sxy)<-sxy
      # proj4string(sxy)<-proj4string(habPolygon)
      # sxy.sf<-st_as_sf(sxy,"POINTS")
      # sxy.sf<-st_intersects(sxy.sf,habPolygon.sf)
      # #sxy<-gIntersection(sxy,habPolygon)
      # length(sxy.sf)
      # 
   # sxy<-posterior.sxy[1,,]
   # out<-apply(posterior.sxy[1:10,,],1,function(sxy){
   #    sxy<-data.frame(sxy)
   #    coordinates(sxy)<-sxy
   #    proj4string(sxy)<-proj4string(habPolygon)
   #    sxy<-gIntersection(sxy,habPolygon)
   #    sxy#length(sxy)
   # })

   
   if(is.null(habRaster) & is.null(habResolution)){
   habPolygon.sf<-st_as_sf(habPolygon,"POLYGON")
   
   i<-1
   counts<-unlist(lapply(1:dim(posterior.sxy)[1],function(i){
      print(i)
      
      sxy<-data.frame(posterior.sxy[i,posterior.z[i,]%in%alive.states,])

      coordinates(sxy)<-sxy
      proj4string(sxy)<-proj4string(habPolygon)
      sxy.sf<-st_as_sf(sxy,"POINTS")
      sxy.sf<-st_intersects(sxy.sf,habPolygon.sf)
      #sxy<-gIntersection(sxy,habPolygon)
      length(sxy.sf)
   }))
   
   out<-c(mean=mean(counts),LCI=as.numeric(quantile(counts,0.025)),UCI=as.numeric(quantile(counts,0.975)))
   return(out)
   }
   #dim(posterior.z)
   
   
   # JOE'S  FANTASTIC BULLSHIT
   # Filter out posterior sxy for dead individuals
   
#SLOWER: deadId <- which(posterior.z != alive.states, arr.ind = TRUE)
   deadId <- which(!posterior.z%in%alive.states, arr.ind = TRUE)

   
   # for(i in 1:dim(deadId)[1]){
   #    posterior.sxy[deadId[i,1],deadId[i,1],1:2] <-
   # }
   post.x <- posterior.sxy[ , ,1]
   post.x[deadId] <- -999
   posterior.sxy[,,1] <- post.x
   
   # Convert habRaster to SpatialGrid object 
   spRaster <- as(habRaster, "SpatialGridDataFrame")

   # Calculate cell-specific density for each iteration
   #curCoords<-posterior.sxy[1,,]
   Density <- apply(X = posterior.sxy, FUN = function(curCoords, inRaster){
      #print("1")
      # Get the cell index for each coordinate
      curGridIndeces <- getGridIndex(curCoords, getGridTopology(inRaster), all.inside = FALSE)
      # Get rid of coordinate indeces that do not fall in a cell
      curGridIndeces <- curGridIndeces[!is.na(curGridIndeces)]
      # Create a frequency table of grid indeces
      freqTable <- table(as.character(curGridIndeces))
      # Initialise an output vector
      outVector <- rep(0, nrow(coordinates(inRaster)))
      outVector[as.integer(names(freqTable))] <- freqTable
      outVector[is.na(spRaster@data$layer)] <- NA
      outVector
   }, inRaster = spRaster, MARGIN = 1)
#summary(outVector)
   
   # median(1:4)
   # plot(density(counts))
   # quantile(c(1,2,3,4),0.025)
   # quantile(sample(4),0.025)
   # quantile(sample(4,5,replace=TRUE),0.025)
   # 

   counts<-apply(Density,2,sum,na.rm=TRUE)
   #quantile(counts,0.025,type=1)#,na.rm=TRUE)
   if(!return.all)out<-c(mean=mean(counts),median=median(counts),mode=GetMode(counts),LCI=as.numeric(quantile(counts,0.025,na.rm=TRUE)),UCI=as.numeric(quantile(counts,0.975,na.rm=TRUE)))
 
   if(return.all)out<-counts
   
   return(out)  
}
   
   
#    # EXTRACT SUMMARY STATISTICS PER HABITAT CELL
#    meanDensity <- apply(Density, 1, mean)
#    CI <- apply(Density, 1, function(x)quantile(x, c(0.025,0.5,0.975), na.rm = TRUE))
#    
#    # FEED IN RASTERS
#    meanDensity.r <- lowCIDensity.r <- upCIDensity.r <- medianDensity.r <- habRaster
#    meanDensity.r[] <- meanDensity
#    lowCIDensity.r[] <- CI[1, ]
#    medianDensity.r[] <- CI[2, ]
#    upCIDensity.r[] <- CI[3, ]
#    
#    if(plot.check){
#       par(mfrow = c(1,2))
#       ## Plot mean density 
#       max.cut <- max(meanDensity.r[], na.rm = TRUE)
#       cuts <- seq(0, max.cut, length.out = 101) 
#       pal <- rev(terrain.colors(length(cuts)))
#       plot(meanDensity.r, breaks = cuts, col = pal, main ="Mean density", legend=FALSE, axes=FALSE)
#       plot(meanDensity.r, legend.only=TRUE, col=pal,
#            legend.width=1, legend.shrink=1,
#            axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
#                           labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
#                           cex.axis=0.6),
#            legend.args = list(text=paste("Density (ind.",(habResolution/1000)^2,"km-2)",sep=""), side=4, font=2, line=2.5, cex=0.8))
#       plot(habPolygon, add = TRUE)
#       
#       ## Plot CI width
#       max.cut <- max(upCIDensity.r[]-lowCIDensity.r[], na.rm = TRUE)
#       cuts <- seq(0, max.cut, length.out = 101) 
#       pal <- rev(heat.colors(length(cuts)))
#       plot(upCIDensity.r-lowCIDensity.r, breaks = cuts, col = pal, main ="95% CI width", legend = FALSE, axes=FALSE)
#       plot(upCIDensity.r, legend.only=TRUE, col=pal,
#            legend.width=1, legend.shrink=1,
#            axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
#                           labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
#                           cex.axis=0.6))
#       plot(habPolygon, add = TRUE)
#    }
#    return(list("mean.Density" = meanDensity.r,
#                "median.Density" = medianDensity.r,
#                "upperCI.Density" = upCIDensity.r,
#                "lowerCI.Density" = lowCIDensity.r))
# }
