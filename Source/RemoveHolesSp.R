#' @title RemoveHolesSp() removes holes from a spatialPolygons.
#'
#' @description
#' \code{RemoveHolesSp} returns a \code{SpatialPolygons} object with no holes
#' 
#' @param sp A \code{sp} SpatialPolygonsDataFrame or a SpatialPolygons
#' @return A \code{SpatialPolygons} 
#' 
#' @examples
#' spNoHoles <- RemoveHolesSp(sp)


RemoveHolesSp <- function(sp){
  
  sp <- as(sp,"SpatialPolygons")

    for(i in 1:length(sp)){
        is.hole <- unlist(lapply( sp@polygons[[i]]@Polygons, function(x) x@hole))
        sp@polygons[[i]]@Polygons <- sp@polygons[[i]]@Polygons[!is.hole]
     }
  return(sp)
  
}