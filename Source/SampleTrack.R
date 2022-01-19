#' @title Function to sample points along a line (track) with constant distance 
#' @description
#' \code{UTMToGrid} returns a dataframe with x and y coordinates.
#' 
#' @param x A \code{SpatialLinesDataFrame} sp class object
#' @param d A \code{numeric}  Sample distance. For regular sample. 
#' @param P A \code{numeric} Proportional sample size (length * p), expected value is 0-1. For regular or random.
#' @param n A \code{numeric} Fixed sample size. For regular or random
#' @param type	A \code{string} Defines sample type. Options are "regular" or "random". A regular sample results in a systematic, evenly spaced sample.
#' @param longlat	\code{logical} TRUE/FALSE is data in geographic units, if TRUE distance is in kilometres
#' @param poly	\code{SpatialPolygons} to clip the tracks (studyArea), otherwise= NULL
#' @param min.samp A \code{numeric}	Minimal number of sample points for a given line (default is 1 point)
#' @keywords SCR data prep
#'
#' @return A \code{SpatialPointsDataFrame}, object containing the sampled points

#' @examples
#' 
#'my.lines <- SampleTrack(tracks, d = 10000, type = "regular", longlat = TRUE, min.samp = 1)
#'
#' 
### this is the function modified from the package spatialEco


SampleTrack <- function(   x
                         , d = 10000
                         , p = NULL
                         , n = NULL
                         , type = "regular"
                         , longlat = FALSE
                         , poly = NULL
                         , min.samp = 1
                         , ...){

   ### this is the function modified from the package spatialEco
   
   if(!is.null(poly)){
      x <- intersect(x, poly)
   }
   
   
  samp.size <- function(l, p.s = NULL, n.s = NULL, d.s = d, longlat.s = FALSE, 
                        min.samp.s = 1) {
    line.length <- sp::SpatialLinesLengths(l, longlat = longlat.s)
    if (!is.null(n.s)) {
      ns = n.s
    }
    else if (!is.null(p.s)) {
      ns <- round(line.length[1] * p.s, digits = 0)
    }
    else {
      ns <- round((line.length[1]/d.s), digits = 0)
    }
    if (ns < 1) 
      ns = min.samp
    return(ns)
  }   
  
store <-list() # store each points in a list to make the function faster 


  if (!inherits(x, "SpatialLinesDataFrame")) 
    stop("x is not a SpatialLinesDataFrame object")
  if ((!is.null(p) == TRUE) && (!is.null(n) == TRUE)) 
    stop("Cannot have both a fixed and proportional sample")
  lids <- rownames(x@data)

  lsub <- x[rownames(x@data) == lids[1], ]
  ns <- samp.size(lsub ,d.s = d)
  lsamp <- sp::spsample(lsub, n = ns, type = type)
  
  if(is.null(lsamp)){
    ns <- 20 ### here is to trick it so when sp sample doesnt manage to sample one single point. putting 20 points help having one point at least 
    lsamp <- sp::spsample(lsub, n = ns, type = type)[1,]
    
    }
  
  ## store in list and then do.call(rbind at the end)
  store[[1]] <- sp::SpatialPointsDataFrame(lsamp, data = data.frame(LID = rep(as.numeric(lids[1]), 
                                                                     length(lsamp))))
  
  
  if (length(lids) > 1) {
    for (i in 2:length(x)) {
      lsub <- x[rownames(x@data) == lids[i], ]
      ns <- samp.size(lsub,d.s = d)
  
      lsamp <- sp::spsample(lsub, n = ns, type = type)
      if(is.null(lsamp)){
        ns <- 20 ### here is to trick it so when sp sample doesnt manage to sample one single point. putting 20 points help having one point at least 
        lsamp <- sp::spsample(lsub, n = ns, type = type)[1,]
        }
      lsamp <- sp::SpatialPointsDataFrame(lsamp, data = data.frame(LID = rep(as.numeric(lids[i]), 
                                                                             length(lsamp))))
      
      store[[i]] <- lsamp
      #print(i)
    }
    
  }
  results <- do.call(rbind, store)## do.call the list to improve perfomance 
  return(results)
}
