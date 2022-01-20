#' @title Function to simulate a spatially autocorelated raster 
#'
#' @description
#' \code{SimulateAutocorrelatedRaster} returns a raster with values being spatially autocorelated.
#' it is based on a multivariate normal distribution with the function being taken from http://rstudio-pubs-static.s3.amazonaws.com/9688_a49c681fab974bbca889e3eae9fbb837.html
#' 
#' 
#' @param sp object containing ths sp file (e.g. the habitat or detector grid) .
#' @param NRaster A \code{numeric} for the number of raster files to simulate. 
#' @param phi A \code{numeric} parameter describing how rapidly the correlation declines with distance.
#' @param scaled A \code{logical} if raster values should be scaled or not 
#' @param focal if an average focal moving window, a \code{numeric} with the size of the moving window to be used. if focal = 5, a 5*5 moving window is used.  
#' @param seed if results to be reproduced, use a \code{numeric} value for the seed.
#' @param plot Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated during simulations.
#' 

#'
#' @examples
#' # Generate a spatial point data frame with simulated activity center locations
#' MyDensityRaster <- SimulateAutocorrelatedRaster(sp=myHabitat.list$habitat.sp
#'                                                  , NRaster=1
#'                                                  , phi = 0.1
#'                                                  , scaled = T
#'                                                  , focal = NULL
#'                                                  , plot = TRUE
#'                                                  , seed=NULL)

SimulateAutocorrelatedRaster_EM <- function(sp
                                           , NRaster = 1
                                           , phi = 1 
                                           , scaled = FALSE
                                           , focal = NULL
                                           , plot = TRUE
                                           , seed = NULL
                                           
                                           
                                           
){
  
  
  #function taken from http://rstudio-pubs-static.s3.amazonaws.com/9688_a49c681fab974bbca889e3eae9fbb837.html
  rmvn <- function(n, mu = 0, V = matrix(1),seed=NULL) {
    p <- length(mu)
    if (any(is.na(match(dim(V), p))))
      stop("Dimension problem!")
    #D <- chol(V)
    D <- V
    if(!is.null(seed)){set.seed(seed)}
    t(matrix(rnorm(n * p, mean=0, sd=1), ncol = p) %*% D + rep(mu, rep(n, p)))
  }
  
  # GET DITANCE MATRIX FOR THE HABITAT
  
  distance <- gDistance(sp, sp, byid = T)
  
  r <- list()
  moran <- 0
  if(length(phi)==1){ phi<- rep(phi,NRaster) } 
  for( i in 1:NRaster){
    # SIMULATE THE COVARIATE 1
    distance1 <- exp(-phi[i] * distance)
    n <- dim(distance1)[1]
    if(!is.null(seed)){seed1 <- seed+i}else{
      seed1 <- seed
    }
    X <- rmvn(   n=1
                 , mu = rep(0, n)
                 , V = exp(-phi[i] * distance)
                 , seed = seed1
    )
    
    r[[i]] <- rasterFromXYZ(cbind(coordinates(sp), X))
    
    if(scaled) r[[i]][] <- scale(r[[i]][])
    
    if(!is.null(focal)){
      ## apply the focal 
      f.rast <- function(x) mean(x, na.rm=T) 
      MovingWindowSize <-  matrix(1, focal, focal)
      r[[i]] <- focal(r[[i]], MovingWindowSize, f.rast, pad=T)
    }
    
    moran[i] <- Moran(r[[i]])
  }
  
  r <- stack(r)
  
  names(r) <- paste("Moran's Index", round(moran, digits = 2))
  
  if(plot){
    plot(r)
  }
  
  return(r)
  
  
  
}
