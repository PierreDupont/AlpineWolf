EstimateN <- function( r.origin                 # A raster file of the study area 
                       , IDCells.mx               # A matrix with the cells ID of the study area  
                       , posterior.sy             # The posterior x coordinates of individual ACs 
                       , posterior.sx             # The posterior y coordinates of individual ACs
                       , posterior.z              # The posterior of z 
                       , alive.states             # Alive states in z
                       , poly                     # Polygon of the study area
                       , time.interval = NULL     # Years to retained 
                       , CI = c(0.025,0.975))     # Confidence INterval to be returned            
{
   ## Set the parameters of the results to be returned
   n.years <- ifelse(length(dim(posterior.z))==2,1,ifelse(is.null(time.interval),dim(posterior.z)[3],length(time.interval)))
   if(is.null(time.interval)){time.interval <- 1:n.years} 
   n.iter <- dim(posterior.z)[1]
   n.individuals <- dim(posterior.z)[2]
   
   ## Filter the posteriors
   sx <- sy <- Z <- array(0,dim = c(n.iter,n.individuals,n.years))
   if(length(dim(posterior.z))!=2){posterior.z <- posterior.z[ , ,time.interval]}
   if(length(dim(posterior.sy))!=2)
      {sy <- posterior.sy[ , ,time.interval]
      sx <- posterior.sx[ , ,time.interval]}
   if(length(dim(posterior.sy))==2)
      {
      for(t in 1:n.years)
         {
         sy[ , ,t] <- posterior.sy
         sx[ , ,t] <- posterior.sx
         }
      }
   
   ## Assign a 1 to all alive individuals
   Z[which(posterior.z %in% alive.states)] <- 1
   
   # Transform AC locations to cell ID
   sx <- trunc(sx)+1
   sy <- trunc(sy)+1
   
   ## Keep the AC locations, Iterations and time of alive individuals only (z == 1)
   sy <- sy[which(Z==1)]
   sx <- sx[which(Z==1)]
   ITER <- which(Z==1,arr.ind=TRUE)[,1] 
   TIME <- which(Z==1,arr.ind=TRUE)[,3] 
   
   ## Identify the cells in the Study Area
   cell.ID.clip <-  raster::values(r.origin)
   if (!is.null(poly))
   {
      raster::values(r.origin) <- IDCells.mx                  # Populate the raster with the ID of the cells. 
   r.clip  <- mask(r.origin, poly)                 # Clip the raster with the polygon.clip 
   val <- raster::values(r.clip)                           # Extract only the cells ID that are inside the polygon clip 
   cell.ID.clip <- val[!is.na(val)]                # Remove NAs
   }
   ##----- Create the matrix of population sizes for each iteration -----
   N.counts <- matrix(0, n.iter, n.years)
   for (x in 1:length(sy))
      {
      if(IDCells.mx[sy[x],sx[x]] %in% cell.ID.clip)
         {N.counts[ITER[x],TIME[x]] <- N.counts[ITER[x],TIME[x]] + 1}
      }
      
   ##----- Produce output list -----
   MEDIAN <- apply(N.counts,2,median)
   MEAN <- apply(N.counts,2,mean)
   CONF.INT.min <- apply(N.counts,2,FUN= function(x){quantile(x,CI[1])})
   CONF.INT.max <- apply(N.counts,2,FUN= function(x){quantile(x,CI[2])})
   AREA <- gArea(poly)
   SUMMARY <- cbind.data.frame( time.interval
                              , MEDIAN
                              , MEAN
                              , CONF.INT.min
                              , CONF.INT.max
                              , rep(AREA,n.years))
   
   names(SUMMARY) <- c("Year", "Median", "Mean", paste(CI[1],"%"), paste(CI[2],"%"),"Area")
   
   return(list(N.counts=N.counts, summary=SUMMARY))
}
