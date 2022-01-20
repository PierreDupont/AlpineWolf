#' @title Function to make a  extract covariates from a rasterstack
#' #'
#' @description
#' \code{ExtractCovariates} returns a \code{matrix}  object with with the length of the spatialpoint object in row and of the raster stack in columns
#' 
#' @param rasterstack A \code{rasterstack} object with covariates to be extracted
#' @param detectors.sp A \code{SpatialPoints} with le location of the covariates to be extracted.
#'                  (ngb method is used, see ?resample from raster). if NUll no resampling is performed and original resolution is retained.
#' @buffer A \code{numeric} with the length of the buffer used to extract the values. if buffer is used, the mean of all cells computed  
#' @examples
#' 
#'  covariates <- ExtractCovariates(  rasterstack = rasterstack
#'                                   , detectors.sp = detectors.sp
#'                                   , buffer= NULL)
#'                                   
#'
#'
#'
#' 
#' 



ExtractCovariates <- function(  rasterstack = rasterstack
                              , detectors.sp = detectors.sp
                              , buffer= NULL){
   
   extval <- raster::extract(rasterstack, detectors.sp, buffer=buffer, fun=mean)
   
   return(extval)
   
}



