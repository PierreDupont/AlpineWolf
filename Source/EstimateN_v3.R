EstimateN_v3 <- function(habRaster = NULL, # RASTER  FOR WHICH REGION ESTIMATES SHOULD BE ESTIMATED. 
                         posterior.sxy,
                         posterior.z,
                         alive.states,
                         return.all = FALSE,
                         regionEstimates = FALSE # IF TRUE= CELL BASED STATISTICS ARE REPORTED. IF FALSE, REGION SPECIFIC ESTIMATES ARE PROVIDED.
){
   
  # IDENTIFY ITERATIONS FOR WHICH INDIVIDUALS ARE NOT ALIVE
   # deadId <- which(!posterior.z%in%alive.states, arr.ind = TRUE)
   # post.x <- posterior.sxy[ , ,1]
   # post.x[deadId] <- -999
   # posterior.sxy[,,1] <- post.x
   
   #---[RB] fixed 2019-08-29
   deadId <- which(!posterior.z %in% alive.states)#, arr.ind = TRUE)
   posterior.sxy[,,1][deadId]<--999
   #hist(apply(posterior.sxy[,,1],1,function(x)sum(x!=-999)))
   
   
   # Convert habRaster to SpatialGrid object 
   spRaster <- as(habRaster, "SpatialGridDataFrame")
   
   # Calculate cell-specific density for each iteration
   Density <- apply(X = posterior.sxy, FUN = function(curCoords, inRaster){
      #print("1")
      # Get the cell index for each coordinate
      curCoords <- curCoords[curCoords[,1]!=-999,]
      curGridIndeces <- getGridIndex(curCoords, getGridTopology(inRaster), all.inside = FALSE)
      # Get rid of coordinate indeces that do not fall in a cell
      curGridIndeces <- curGridIndeces[!is.na(curGridIndeces)]
      # Create a frequency table of grid indeces
      freqTable <- table(as.character(curGridIndeces))
      # Initialise an output vector
      outVector <- rep(0, length(inRaster))
      outVector[as.integer(names(freqTable))] <- freqTable
      outVector
   }, inRaster = spRaster, MARGIN = 1)
   #RETURN NA IF OUTSIDE OF HABITAT 
   Density[is.na(spRaster@data$layer),] <- NA
   
   ## IF REGION ESTIMATES ARE ASKED  
   if(regionEstimates){
      PosteriorsRegion <- list()
      IDregions <- unique(spRaster@data$layer) 
      IDregions <- IDregions[!is.na(IDregions)]
      IDregions <- as.character(IDregions)
      
      summary <- matrix(NA, nrow=length(IDregions)+1, ncol=5)
      rownames(summary) <- c(IDregions, "Total")
      colnames(summary) <- c("mean","median","mode", "95%CILow","95%CIHigh")
      
      for(i in IDregions){
         temp <- Density
         temp[!(spRaster@data$layer %in% i),] <- NA
         PosteriorsRegion[[i]] <- apply(temp,2,sum,na.rm=TRUE)
         summary[which(row.names(summary)==i),] <- c(round(mean(PosteriorsRegion[[i]]),digits=2),
                                                     median(PosteriorsRegion[[i]]),
                                                     GetMode(PosteriorsRegion[[i]]),
                                                     quantile(PosteriorsRegion[[i]], probs=c(0.025,0.975)))
      }
      
     
      
      # GET TOTAL ESTIMATES 
      TotalPosteriorsRegion <- apply(Density,2,sum,na.rm=TRUE)
      PosteriorsRegion[["ALL"]] <-  TotalPosteriorsRegion
      summary[length(IDregions)+1,] <- c(round(mean(TotalPosteriorsRegion),digits=2),
                                         median(TotalPosteriorsRegion),
                                         GetMode(TotalPosteriorsRegion),
                                         quantile(TotalPosteriorsRegion, probs=c(0.025,0.975)))
      
      out <- list(PosteriorsRegion = PosteriorsRegion,
                  summary = summary)
      ## IF CELL POSTERIORS ARE ASKED  
   }else{
      out <- list(
      PosteriorsCellsMean = apply(Density, 1, mean),
      PosteriorsCellsMedian = apply(Density, 1, median),
      PosteriorsCellsCILow = apply(Density, 1, function(x) quantile(x, c(0.025), na.rm = TRUE)),
      PosteriorsCellsCILow = apply(Density, 1, function(x) quantile(x, c(0.975), na.rm = TRUE))
      )
     
    return(out)

   }
   
}

