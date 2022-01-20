#' @title Function to create a transition layer for least cost path analysis.
#'
#' @description
#' \code{SimulateDetection} returns a list object with \code{} 
#' 
#' @param r Raster representing a digital elevation model (or other image that can be interpreted as a cost surface).
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/MakeTransitionLayer
#' @keywords simul
#'
#' @examples
#' # ...
#' ...
#' 
#' 
#' 

MakeTransitionLayer<-function(r,plot=TRUE){
   
   heightDiff <- function(x){x[2] - x[1]}
    
   hd <- transition(r,heightDiff,8,symm=FALSE)
   slope <- geoCorrection(hd, scl=FALSE)
   
   adj <- adjacent(r, cells=1:ncell(r), pairs=TRUE, directions=8)
   speed <- slope
   speed[adj] <- exp(-3.5 * abs(slope[adj] + 0.05))#--For now: use basic cost parameter ("hiker")
   if(plot){
      #-----HILLSHADE
      slope <- terrain(dem.area, opt='slope')
      aspect <- terrain(dem.area, opt='aspect')
      hill <- hillShade(slope, aspect, 65, 270)
      plot(hill, col=grey(0:100/100))
      image(r,col=terrain.colors(25, alpha=0.35), add=TRUE)
   }
   x <- geoCorrection(speed, scl=FALSE)
   return(x)
   
}