#' Spatial Group detection process 
#'
#' The \code{dgroupDet} distribution is a NIMBLE custom distribution which can be
#' used to model and simulate observations (x) of multiple individuals from the 
#' same group over a set of traps defined by their coordinates \emph{trapCoords}.
#' The distribution assumes that all individuals detected together belong to the
#' same group. 
#' 
#' 
#' @name dgroupDet
#'
#' @param x Vector of trap-specific detection frequencies. 
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param size Vector of the number of trials (zero or more) for each trap (\emph{trapCoords}).
#' @param p Baseline detection probability (scalar) used in the half-normal detection function. For trap-specific baseline detection probabilities use argument \emph{p0Traps} (vector) instead.
#' @param alpha Number of traps with at least one detection recorded in \emph{x}; from the \emph{detNums} object returned by the \code{\link{getSparseY}} function. 
#' @param indicator Binary argument specifying whether each group is available for detection (indicator = 1) or not (indicator = 0).
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#'
#' @return The log-likelihood value associated with the vector of detections, given the location of the activity center (s),
#'  and the half-normal detection function : \eqn{p = p0 * exp(-d^2 / 2 \sigma^2)}.
#'
#' @author Pierre Dupont
#'
#' @import nimble
#' @importFrom stats dbinom rbinom pbinom 
#'
#' @examples
#' ## Create habitat grid cell coordinates
#' xBounds <- c(0.5,25.5)
#' yBounds <- c(0.5,50.5)
#' habResolution <- 1
#' xCoords <- seq(xBounds[1],xBounds[2],habResolution)
#' yCoords <- seq(yBounds[1],yBounds[2],habResolution)
#' coordsHabitatGridCenter <- expand.grid(
#'   list(x = xCoords,
#'          y = yCoords))
#'          colnames(coordsHabitatGridCenter) <- c("x","y")
#'          numHabWindows <- nrow(coordsHabitatGridCenter)
#'          
#'          ## Create trap coordinates
#'          bufferWidth <- 2
#'          trapResolution <- 0.7
#'          trapCoords <- expand.grid(
#'            list(x = seq(xBounds[1]+bufferWidth,xBounds[2]-bufferWidth,trapResolution),
#'                  y = seq(yBounds[1]+bufferWidth,yBounds[2]-bufferWidth,trapResolution)))
#'                  colnames(trapCoords) <- c("x","y")
#'                  numTraps <- nrow(trapCoords)
#'                  oper <- rep(1,numTraps)
#'                  
#'                  ## RESCALE COORDINATES 
#'                  scaledCoords <- scaleCoordsToHabitatGrid( coordsData = trapCoords,
#'                                                            coordsHabitatGridCenter = coordsHabitatGridCenter)
#'                                                            scaledtrapCoords <- scaledCoords$coordsDataScaled
#'                                                            habitatMask <- matrix(1, nrow = length(yCoords), ncol = length(xCoords), byrow = TRUE)
#'                                                            
#'                                                            
#'  ## CREATE LOCAL OBJECTS 
#'  localTraps <- getLocalObjects( habitatMask = habitatMask,
#'                                coords = scaledtrapCoords,
#'                                dmax = 10,
#'                                resizeFactor = 1,
#'                              plot.check = TRUE)
#' 
#' ## simulate population
#' N <- 50
#' M <- 150
#' indicator <- c(rep(1,N),rep(0,M-N))
#' s <- cbind(runif(n = M, min = xBounds[1], max = xBounds[2]),
#'            runif(n = M, min = yBounds[1], max = yBounds[2]))
#' colnames(s) <- c("x","y")
#' scaledS <- scaleCoordsToHabitatGrid( coordsData = s,
#'                                      coordsHabitatGridCenter = coordsHabitatGridCenter)$coordsDataScaled
#' 
#' ## PLOT CHECK (RESCALE SCALE)
#' plot(scaledCoords$coordsHabitatGridCenterScaled[,"y"]~scaledCoords$coordsHabitatGridCenterScaled[,"x"]) 
#' points(scaledtrapCoords[,"y"]~scaledtrapCoords[,"x"], col = "blue",pch=3) 
#' points(scaledS[,"y"]~scaledS[,"x"],col=adjustcolor("red",0.3),pch=16)
#' points(scaledS[indicator>0,"y"]~scaledS[indicator>0,"x"],col="red",pch=16) 
#' 
#' 
#' ## DETECTION PARAMETERS
#' lambda0 <- 1
#' sigma <- 1.5
#' 
#' # III. USING THE RANDOM GENERATION FUNCTION 
#' y <- rgroupDet( n = 1,
#'                           lambda0 = lambda0,
#'                           sigma = sigma,
#'                           s = scaledS,
#'                           trapCoords = scaledtrapCoords,
#'                           oper = oper,
#' indicator = indicator)
#' 
#' ## Plot number of detections per trap
#' points( scaledtrapCoords[,"y"]~scaledtrapCoords[,"x"],
#'         pch = 21,
#'         bg = adjustcolor("black", 0.5),
#'         col = "black",
#'         cex = y/2)
#' 
#' 
#' # II. USING THE DENSITY FUNCTION 
#' system.time(
#' logProb <-  dgroupDet( x = y,
#'                   lambda0 = lambda0,
#'                   sigma = sigma,
#'                   s = s,
#'                   trapCoords = trapCoords,
#'                   oper = oper,
#'                   indicator = indicator,
#'                   log = T)
#'                    )
#' print(logProb)
#' 
#' @export
NULL

#' @rdname dgroupDet
#' @export
dgroupDet <- nimbleFunction(
  run = function( x = double(0),           ## data: Number of wolves detected together 
                  size = double(1),        ## Group sizes
                  p = double(1),           ## Group-specific detection probabilities
                  alpha = double(0),       ## Cohesion (i.e. prob of detection given visit)
                  indicator = double(1),   ## Group-specific augmentation indicator
                  log = integer(0, default = 0)
  ){
    ##-- Return type
    returnType(double(0))
    
    ##-- Derive probability of 0 visit:
    numGroups <- length(p)
    sumP <- sum(p[1:numGroups])
    pNULL <- exp(-sumP)
    risk <- 1 - pNULL
    
    ##-- if n = 0, return probability of 0 visit :
    if(x == 0){
      if(log == 0){return(pNULL)}else{return(log(pNULL))}
    }
    
    ##-- If x > 0; loop over groups and calculate the probability
    ##-- to detect X individuals simultaneously given the group size 
    ##-- and group-specific probability of visit :
    totalProb <- 0
    for(g in 1:numGroups){
      if(indicator[g] > 0){
        if(x <= size[g]){
          ##-- Probability of visit by group[g]
          probThisGroup <- risk * p[g]/sumP
          
          ##-- Probability of x wolves from group[g] detected together 
          ##-- (Truncated binomial distribution)
          logNormCst <- pbinom( q = 0,
                                size = size[g],
                                prob = alpha,
                                lower.tail = 0,
                                log.p = 1)
          
          logProbThisNum <- dbinom( x = x,
                                    size = size[g],
                                    prob = alpha,
                                    log = 1) -  logNormCst
          
          ##-- Probability of the data
          totalProb <- totalProb + probThisGroup * exp(logProbThisNum)
        }#if group size is less than the number of observed individuals
      }#if group is not augmented
    }#g
    if(log) return(log(totalProb)) else return(totalProb)
  })

    



#' @rdname dgroupDet
#' @export
rgroupDet <- nimbleFunction(
  run = function( n = double(0, default = 1),   
                  size = double(1),        
                  p = double(1),         
                  alpha = double(0),       
                  indicator = double(1)    
  ){
    ##-- Return type
    returnType(double(0))
    
    ##-- Derive probability of 0 visit:
    numGroups <- length(p)
    sumP <- sum(p[1:numGroups])
    pNULL<- exp(-sumP)
    risk <- 1 - pNULL
    
    ##-- if n = 0, return 0
    numVisits <- rbinom(n = 1, size = 1, prob = risk)
    if(numVisits == 0){return(numVisits)}
    
    ##-- If x > 0; loop over groups and calculate the probability
    ##-- to detect X individuals simultaneously given the group size 
    ##-- and group-specific probability of visit 
    whichGroup <- rcat(n = 1, prob = p[1:numGroups]*indicator[1:numGroups])
    
    out <- rbinom( n = 1,
                   size = size[whichGroup],
                   prob = alpha)
    
    return(out)
  })


## dbinomLocal_normal_SEM
registerDistributions(
  list(
    dgroupDet = list(
      BUGSdist ='dgroupDet(size, p, alpha, indicator)',
      Rdist = c('dgroupDet(size, p, alpha, indicator)'),
      types = c('value = double(0)',
                'size = double(1)', 'p = double(1)',
                'alpha = double(0)', 'indicator = double(1)'),
      discrete = TRUE,
      mixedSizes = TRUE,
      pqAvail = FALSE)
  ), verbose = F)

