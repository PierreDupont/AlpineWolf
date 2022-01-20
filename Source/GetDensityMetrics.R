#' @title Get comparisons true-estimated density metrics
#'
#' @description
#' \code{GetDensityMetrics} returns a \code{list} object with diferent metrics about the comparisons of estimated spatial and estimated density
#'
#' @param AC.sp A \code{SpatialPoints} object containig the locations of the ACS from simulated individuals 
#'	@param myHabitat.list A \code{list} object from makeHabitat function 
#' @param true.density.sp A \code{raster} object from SimulateTrueDensity() function that estimates true density. 
#' @param predicted.density.sp A \code{raster} object from PlotPredictedACs containing the estimated density from jagsmodel 
#' @param jagsoutput A \code{list} object from the output of the model.
#' @param detected A \code{detected} A vector with the ID of individuals that were detected. this ID should match with the ID of the AC.sp. (this is used to calculate distance between TRUE and Estimated ACS)
#' @param methods.correlation A \code{character} indicating which correlation coefficient is to be used. One of "pearson", "kendall", or "spearman"
#' @param ngb neighborhood size. Either a single integer or a vector of two integers c(nrow, ncol), see corLocal for details 
#' @return A \code{list} object metrics about the performance of the density estimates.
#' $distance.ACS:  return metrics about distances from true to estimated ACS. 
#'                   $summary         = metrics summarizing distance between estimated and predicted ACs for all detected individuals 
#'                   $individual.dist = metrics summarizing distance between estimated and predicted ACs for each detected individuals 
#' $correlation : returns correlations between predicted and true density
#'                local.correlation = returns Local correlation coefficient for two RasterLayer objects 
#'                                    (using a focal neighborhood) see corLocal for details 
#'                 cor.test = returns output form cor.test between the two raster layers                   
#' $bias: calculates bias metrics between predicted.density.sp and true.density.sp
#'           $density.bias.sp = a \code{raster} object with bias between estimated and predicted bias. 
#'           $mean.bias       = mean bias 
#'           $sd.bias         = standard deviation of the bias. 
#'           $RMSE            = return the root mean square error
#'           $NRMSE           = return the normalized root mean square error 
#'                              $maxmin = using maxmin methods  RMSE/max(density.bias.vec, na.rm=T) -   min(density.bias.vec, na.rm=T)
#'                              $sd     = using sd methods      RMSE/sd(density.bias.vec, na.rm=T)
GetDensityMetrics <- function( Pred.Dens
                             , Sim.Dens
                             , AC.sp
                             , Detectors = NULL
                             , y = NULL
                             , plot.check = TRUE
                             , poly = NULL
                             , cookie.cut = FALSE)
   {
   pred.dens <- Pred.Dens$mean
   sim.dens <- Sim.Dens$POP.r
   
   ## ==== 1.COOKIE CUT TO THE PROVIDED POLYGON ====
   if(cookie.cut){
      if(is.null(poly)){stop("THE POLYGON TO SUBSET THE DENSITY RASTER IS MISSING!")}
      pred.dens <- crop(pred.dens, poly)
      sim.dens <- crop(sim.dens, poly)
      }
   
   ## ==== 2.CALCULATE DIFERENCES & BIASES ====
   diff.dens <- pred.dens - sim.dens
   mean.diff.dens <- mean(diff.dens[],na.rm=TRUE)
   sd.diff.dens <- sd(diff.dens[],na.rm=TRUE)
   
   bias.dens <- abs(diff.dens)
   mean.bias.dens <- mean(bias.dens[],na.rm=TRUE)
   sd.bias.dens <- sd(bias.dens[],na.rm=TRUE)
   
   relative.diff.dens <- diff.dens/sim.dens
   mean.relative.diff.dens <- mean(relative.diff.dens[],na.rm=TRUE)
   sd.relative.diff.dens <- sd(relative.diff.dens[],na.rm=TRUE)
   
   relative.bias.dens <- abs(relative.diff.dens)
   mean.relative.bias.dens <- mean(relative.bias.dens[],na.rm=TRUE)
   sd.relative.bias.dens <- sd(relative.bias.dens[],na.rm=TRUE)
   
   ## ==== 3.PLOTS ====
   if(plot.check == TRUE){
      par(mfrow = c(1,3))
      alldens <- range(na.omit(c(pred.dens[],sim.dens[])))
      max.cut <- max(alldens, na.rm = TRUE)
      cuts <- seq(0, max.cut, length.out = 100) 
      pal <- rev(terrain.colors(length(cuts)))
      
      ## Plot simulated density & all individuals
      plot(sim.dens, breaks = cuts, col = pal, main = "Simulated density")
      if(!is.null(poly))plot(poly, add = TRUE)
      points(AC.sp, pch=16, col="red")
      
      ## Plot estimated density & detected individuals
      plot(pred.dens, breaks = cuts, col = pal, main = "Predicted density")
      if(!is.null(poly))plot(poly, add = TRUE)
      
      if(!is.null(y) & is.null(Detectors)){
         ID.detected <- which(apply(y, 1, function(x)any(x > 0)))
         points(AC.sp[ID.detected, ], col = "blue", pch = 16)
         }
      
      if(!is.null(y) & !is.null(Detectors)){
         ID.detected <- which(apply(y, 1, function(x)any(x > 0)))

        for(i in ID.detected){
           tmp <- Detectors[which(y[ID.detected[i], ] > 0, arr.ind = TRUE), ]
            
            if(length(tmp)==1){                                        # if 1 detection = 1 point only
               points(tmp, col = "blue", pch=16)# plot the point
               }
               
            if(length(tmp)>1){                                        # if more than one sample = plot the ""spider" 
               co <- coordinates(tmp)
               Me <- colMeans(co)                                     # get the mean locations for the ACS
               AC <- SpatialPoints(cbind(Me[1], Me[2]))
               points(AC, col = "blue", pch=16) # plot the estimated ACs
               }#if
            }#i
         }#if
      
      ## Plot spatially-explicit density bias
      plot(bias.dens,  main = "Density difference", col = rev(heat.colors(10)))   
      if(!is.null(poly))plot(poly, add = TRUE)
     
   }
   
   ## ==== 4.CALCULATE RMSE METRICS BETWEEN ESTIMATED AND PREDICTED ACS ==== 
   ## Retrieve density bias summary statistics 
   density.bias.vec <- as.vector(bias.dens)
   # mean
   mean.bias <- mean(density.bias.vec, na.rm = TRUE)
   # standard deviation
   sd.bias <- sd(density.bias.vec, na.rm = TRUE)
   
   ## ==== 5.OUTPUT LIST ==== 
   return(list( density.bias.sp = bias.dens
              , mean.bias = mean.bias
              , sd.bias = sd.bias))
   }


