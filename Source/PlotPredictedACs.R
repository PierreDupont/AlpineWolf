#' @title PlotPredictedACs
#' 
#' @description
#' \code{PlotPredictedACs} Plot individual Predicted ACs locations, detections and locations of simulated ACs
#'
#'
#'
#' @param  matrix.input A matrix of the study area. it can be any matrix as long as it has the same dimensions than the raster used in the analysis
#' @param  raster.origin A raster file with the original projection of the same dimension than matrix input 
#' @param posterior.sy A matrix of posterior sy of  jags. is typically a matrix with individuals (i) in columns i.e. sy[,i] 
#' @param posterior.sx A matrix of posterior sx of jags. is typically a matrix with individuals (i) in columns i.e. sx[,i] 
#' @param posterior.z A matrix of posteriors z of  jags. is typically a matrix with individuals (i) in columns i.e. z[,i] 
#' @param y A matrix of detections, can be a list of matrices when multiple detections from different sources (different y)
#' @param detectors A SpatialPointDataFrame  object with the detectors. can be a list of SpatialPointDataFrame when multiple detections from different sources (different y)
#' @param polygon A SpatialPolygon of the study area 
#' @param IDs A vector with the individual to be plotted. 1 represents first individual in the Y object i.e.  (IDs= c(1:3, 7))
#' @param alive.state a vector of numeric with the alive states  
#' @param Simulated.ACs A SpatialPointsDataFrame containing the simulated ACs positions
#' @param probability A logical Whether maps should be returned as probability of presence (TRUE) or counts (FALSE)
#' @param return.raster.stack A logical whether you also want to return the rasters as a raster stack
#' @param plot.maps  A logical whether you want to plot the graphics or not 
#' @param plot.segments  A logical whether segments and activity centers should be plotted
#' @param plot.ACS  A logical whether (observed) ACs should be plotted or not
#' @param col.detections A vector of color for the locations of the detections
#' @param col.segments A vector of color of the segments linking AC and detection 
#' @param col.ACS A vector of color of the observed ACS 
#' @param sp A SpatialPointsDataFrame object with multiple points (or not) for each individiduals #
#'
#'
#' @return Plot individual Predicted ACs locations, detections and locations of simulated ACs
#' @usage 
#' myPredictedAC.map <- PlotPredictedACs(matrix.input = myHabitat.list$habitat.mx
#'                                       , raster.origin = myHabitat.list$habitat.r
#'                                       , posterior.sy = jagsoutput$JAGS.output$sims.list$sxy[,,2]
#'                                       , posterior.sx = jagsoutput$JAGS.output$sims.list$sxy[,,1]
#'                                       , posterior.z = jagsoutput$JAGS.output$sims.list$z
#'                                       , y = y
#'                                       , detectors = myDetectors$detector.sp
#'                                       , polygon = myStudyArea.poly
#'                                       , IDs = c(1:dim(y)[1])
#'                                       , probability = TRUE 
#'                                       , Simulated.ACs = mySimulatedACs[detected,]
#'                                       , return.raster.stack = TRUE
#'                                       , plot.maps = TRUE
#'                                       , plot.segments = TRUE
#'                                       , plot.ACS = TRUE
#'                                       , col.detections = "black"
#'                                       , col.segments = "blue"
#'                                       , col.ACS = "red")
#' @export
#' @examples



PlotPredictedACs <- function(
  matrix.input = matrix.input,    # A matrix of the study area. it can be any matrix as long as it has the same dimensions than the raster used in the analysis
  raster.origin = raster.origin,  # A raster file with the original projection of the same dimension than matrix input 
  posterior.sy = posterior.sy,    # Posterior sy of  jags. is typically a matrix with individuals (i) in columns i.e. sy[,i] 
  posterior.sx = posterior.sx,    # Posterior sx of jags. is typically a matrix with individuals (i) in columns i.e. sx[,i] 
  posterior.z = posterior.z,      # Posteriors z of  jags. is typically a matrix with individuals (i) in columns i.e. z[,i] 
  y = y,                          # Detection Array
  detectors = detectors,          # A sp object with the detectors 
  polygon = polygon,              # A sp SpatialPolygon of the study area 
  IDs = c(1:3,5,6),               # Individual to be plotted. 1 represents first individual in the Y object. IDs can also be a vector i.e.  (IDs= c(1:3, 7))
  alive.state = c(1),             # Alive states  
  # Additional parameters 
  Simulated.ACs = NULL,            # A SpatialPointsDataFrame containing the simulated ACs positions
  probability =TRUE,              # Whether maps should be returned as probability of presence (TRUE) or counts (FALSE)
  return.raster.stack = TRUE,     # wether you also want to return the rasters as a raster stack
  plot.maps = TRUE,               # whether you want to plot the graphics or not 
  plot.segments = TRUE,           # whether segments and activity centers should be plotted
  plot.ACS = TRUE,                # whether  (observed) ACs should be plotted or not
  col.detections = "black",       # color of the detection points 
  col.segments = "black",         # color of the segments linking AC and detection 
  col.ACS = "red")                # color of the obsvered ACS 
{

 # Replace the alive states by a 1 
   posterior.z[!posterior.z %in% alive.state ]  <- 0
   posterior.z[ posterior.z %in% alive.state ]  <- 1
   

   

  for(i in 1:length(IDs)){       # loop over individuals 
       
      mat <- matrix.input         # Create empty matrix to store the results of the iterations
      mat[] <- 0
    
      # -----  Find the cell based on scaled coordinates ------  
      #j<-15
      for ( j in 1:length(posterior.sy[,IDs[i]])){
      coo <- c(trunc(posterior.sy[j, IDs[i]])+1, trunc(posterior.sx[j,IDs[i]])+1)
      # multiply it by the posterior z to get if individuals were "alive" or not
      mat[coo[1],coo[2]] <-  mat[coo[1], coo[2]] + (1 * posterior.z[j, IDs[i] ])
    }#j
    # obtain probabilities and not counts
    if(probability==TRUE){
      mat  <- mat[]/sum(mat)
    }
    
    # make a raster of the matrix # 
    r.AC <- raster(mat)
    extent(r.AC) <- extent(raster.origin)
    res(r.AC) <- res(raster.origin)
    
    ## ----- Ploting the maps ----- 
    if(plot.maps==TRUE){
      plot(r.AC)
      plot(polygon, add=TRUE)
      # Add the detectors 
      if(!is.list(y)){
      ded <- detectors[which(y[IDs[i],]>0,arr.ind = T),]
    
      # When only one detector, no observed AC
      if(length(ded)==1){
      ded$ID <- 1
      PlotActivityCenters(ded,ID =ded$ID, cex=1, col.detections = col.detections, col.ACS = col.ACS,
                      col.segments = col.segments, plot.segments=TRUE, plot.ACS=TRUE)
    }
    # When more than 1 detector, add segments 
    if(length(ded) >1) {
      ded$ID <- 1
      PlotActivityCenters(ded,ID =ded$ID, cex=1, col.detections = col.detections, col.ACS = col.ACS, col.segments = col.segments,
                      plot.segments=plot.segments, plot.ACS=plot.ACS )
      }
    }else{
       for(s in 1:length(y)){
       ded <- detectors[[s]][which(y[[s]][IDs[i],]>0,arr.ind = T),]
       
       # When only one detector, no observed AC
       if(length(ded)==1){
          ded$ID <- 1
          PlotActivityCenters(ded,ID =ded$ID, cex=1, col.detections = col.detections[s], col.ACS = col.ACS[s],
                              col.segments = col.segments[s], plot.segments=TRUE, plot.ACS=TRUE)
       }
       # When more than 1 detector, add segments 
       if(length(ded) >1) {
          ded$ID <- 1
          PlotActivityCenters(ded,ID =ded$ID, cex=1, col.detections = col.detections[s], col.ACS = col.ACS[s],
                              col.segments = col.segments[s], plot.ACS=plot.ACS )
       }
     }  
    }
    }
       
    # When simulated ACs are provided add them on the plot
    if(i<=length(Simulated.ACs) ){
    if(!is.null(Simulated.ACs)){points(Simulated.ACs[i, ], col="blue", pch=16)}
    }
    ## ----- Return rasterstack ----- 
    if(return.raster.stack == TRUE){
    # stack raster if more than  nb of indiv >1
      if(i>1 ){r.AC1 <- stack(r.AC1, r.AC)}# stack raster of each ID together 
       else{r.AC1 <- r.AC}
     }
  }#i
  ## change names of the raster stack layers, they are numbered from 1 to N 
  if(return.raster.stack == TRUE){
    names(r.AC1) <- as.character(IDs)
    return(r.AC1)
  }
}