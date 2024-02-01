#' Local evaluation of a Poisson Spatial Count detection process 
#'
#' The \code{dpoisLocalSC_normal} distribution is a NIMBLE custom distribution which can be used to model 
#' and simulate Poisson observations (x) of multiple individuals over a set of traps defined by their coordinates \emph{trapCoords}
#' the distribution assumes that an individual’s detection rate at any trap follows a half-normal function of the distance between 
#' the individual's activity center (s) and the trap location. 
#' All coordinates (\code{s} and \code{trapCoords}) should be scaled to the habitat (see \code{\link{scaleCoordsToHabitatGrid}})
#' 
#' The \code{dpoisLocalSC_normal} distribution incorporates two features to increase computation efficiency (see Turek et al., 2021 <doi.org/10.1002/ecs2.3385>  for more details):
#' \enumerate{
#' \item A local evaluation of the individual detection rates calculation (see Milleret et al., 2019 <doi:10.1002/ece3.4751> for more details)
#' \item An indicator (\emph{indicator}) to shortcut calculations for individuals unavailable for detection.
#' }
#' 
#' The \code{dpoisLocalSC_normal} distribution requires x- and y- detector coordinates (\emph{trapCoords}) and activity centers coordinates (\emph{s}) to be scaled to the habitat grid (\emph{habitatGrid}) using the (\code{\link{scaleCoordsToHabitatGrid}} function.)

#' 
#' @name dpoisLocalSC_normal 
#'
#' @param x Vector of trap-specific detection frequencies. 
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param detIndices Vector of indices of traps where the detections in \emph{x} were recorded; from the \emph{detIndices} object returned by the \code{\link{getSparseY}} function. 
#' This argument should not be specified when x is provided as the  \emph{yCombined} object (returned by \code{\link{getSparseY}} ) and when detection data are simulated.
#' @param detNums Number of traps with at least one detection recorded in \emph{x}; from the \emph{detNums} object returned by the \code{\link{getSparseY}} function. 
#' This argument should not be specified when the \emph{yCombined} object (returned by \code{\link{getSparseY}}) is provided as \emph{x} and when detection data are simulated.
#' @param size Vector of the number of trials (zero or more) for each trap (\emph{trapCoords}).
#' @param p0 Baseline detection probability (scalar) used in the half-normal detection function. For trap-specific baseline detection probabilities use argument \emph{p0Traps} (vector) instead.
#' @param p0Traps Vector of baseline detection probabilities for each trap used in the half-normal detection function. When \emph{p0Traps} is used, \emph{p0} should not be provided. 
#' @param sigma Scale parameter of the half-normal detection function.
#' @param s Individual activity center x- and y-coordinates scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).
#' @param trapCoords Matrix of x- and y-coordinates of all traps scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).
#' @param localTrapsIndices Matrix of indices of local traps around each habitat grid cell, as returned by the \code{\link{getLocalObjects}} function.
#' @param localTrapsNum  Vector of numbers of local traps around all habitat grid cells, as returned by the \code{\link{getLocalObjects}} function.
#' @param resizeFactor Aggregation factor used in the \code{\link{getLocalObjects}} function to reduce the number of habitat grid cells to retrieve local traps for.
#' @param habitatGrid Matrix of local habitat grid cell indices, from \emph{habitatGrid} returned by the \code{\link{getLocalObjects}} function. 
#' @param indicator Binary argument specifying whether the individual is available for detection (indicator = 1) or not (indicator = 0).
#' @param lengthYCombined The length of the x argument when the (\emph{yCombined}) format of the detection data is provided;  from the \emph{lengthYCombined} object returned by \code{\link{getSparseY}}
#' 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#'
#' @return The log-likelihood value associated with the vector of detections, given the location of the activity center (s),
#'  and the half-normal detection function : \eqn{p = p0 * exp(-d^2 / 2 \sigma^2)}.
#'
#' @author Pierre Dupont
#'
#' @import nimble
#' @importFrom stats dpois rpois
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
#' ## Create trap coordinates
#' bufferWidth <- 2
#' trapResolution <- 0.7
#' trapCoords <- expand.grid(
#'        list(x = seq(xBounds[1]+bufferWidth,xBounds[2]-bufferWidth,trapResolution),
#'             y = seq(yBounds[1]+bufferWidth,yBounds[2]-bufferWidth,trapResolution)))
#' colnames(trapCoords) <- c("x","y")
#' numTraps <- nrow(trapCoords)
#' oper <- rep(1,numTraps)
#'                  
#' ## RESCALE COORDINATES 
#' scaledCoords <- scaleCoordsToHabitatGrid( coordsData = trapCoords,
#'                                           coordsHabitatGridCenter = coordsHabitatGridCenter)
#' scaledTrapCoords <- scaledCoords$coordsDataScaled
#' habitatMask <- matrix(1, nrow = length(yCoords), ncol = length(xCoords), byrow = TRUE)
#'                                                            
#'  ## CREATE LOCAL OBJECTS 
#'  localTraps <- getLocalObjects( habitatMask = habitatMask,
#'                                coords = scaledTrapCoords,
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
#' y <- rpoisLocalSC_normal( n = 1,
#'                           lambda0 = lambda0,
#'                           sigma = sigma,
#'                           s = scaledS,
#'                           trapCoords = scaledtrapCoords,
#'                           oper = oper,
#'                           localTrapsIndices = localTraps$localIndices,
#'                           localTrapsNum = localTraps$numLocalIndices,
#'                           resizeFactor = 1,
#'                           habitatGrid = localTraps$habitatGrid,
#'                           indicator = indicator)
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
#' # OPTION 1: USING THE RANDOM GENERATION FUNCTIONNALITY 
#' system.time(
#'   logProb <- dpoisLocalSC_normal( x = y,
#'                        lambda0 = lambda0,
#'                        sigma = sigma,
#'                        s = scaledS,
#'                        trapCoords = scaledTrapCoords,
#'                        oper = oper,
#'                        localTrapsIndices = localTraps$localIndices,
#'                        localTrapsNum = localTraps$numLocalIndices,
#'                        resizeFactor = 1,
#'                        habitatGrid = localTraps$habitatGrid,
#'                        indicator = indicator,
#'                       log = T)
#' )
#' print(logProb)
#' 
#' @export
NULL

#' @rdname dpoisLocalSC_normal
#' @export
dpoisLocalSC_normal <- nimbleFunction(
  run = function( x = double(1),
                  lambda0 = double(0),
                  sigma = double(0),
                  s = double(2),
                  trapCoords = double(2),
                  oper = double(1, default = 1.0),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),
                  indicator = double(1),
                  log = integer(0, default = 0)
  ){
    ## Specify return type
    returnType(double(0))
    
    ## Identify dimensions
    numIndividuals <- nimDim(s)[1]
    numTraps <- nimDim(trapCoords)[1]
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma)
    
    ## Calculate trap-specific detection rates
    lambda <- nimNumeric(numTraps)
    ## Loop over individuals 
    for(i in 1:numIndividuals){
      if(indicator[i] > 0){
        ## Retrieve the habitat grid cell ID for this individual's AC 
        sID <- habitatGrid[trunc(s[i,2]/resizeFactor)+1, trunc(s[i,1]/resizeFactor)+1]
        ## Retrieve the local traps surrounding the selected habitat grid cell
        theseLocalTraps <- localTrapsIndices[sID,1:localTrapsNum[sID]]
        ## Loop over local traps only
        for(r in 1:localTrapsNum[sID]){
          ## Calculate distance between AC and local trap
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[i,1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[i,2], 2)
          ## Calculate Poisson rate at this local trap 
          ## This corresponds to the sum of individual lambdas, which follow half-normal distribution
          lambda[theseLocalTraps[r]] <- lambda[theseLocalTraps[r]] + lambda0 * exp(alpha * d2)
        }#r
      }#if
    }#i
    
    ## Calculate the log-probability of the vector of detections
    logProb <- 0.0 
    ## Loop over traps 
    for(j in 1:numTraps){
      logProb <- logProb + dpois( x = x[j],
                                  lambda = lambda[j]*oper[j],
                                  log = TRUE)
    }#J
    # for(j in 1:J) {
    #   ##-- Loop over individuals
    #   for(i in 1:M){
    #     distsq[i,j] <- (s[i,1] - trapCoords[j,1])^2 + (s[i,2] - trapCoords[j,2])^2
    #     lambda[i,j] <- lambda0 * exp(-distsq[g,j] / (2*sigma^2)) * z[g] 
    #   }#i
    #   bigLambda[j] <- sum(lambda[1:M,j]) 
    #   logProb[j] <- dpois(bigLambda[j])
    #   n[j] ~ dpois(bigLambda[j]*oper[j])
    # }#j
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })


#' @rdname dpoisLocalSC_normal
#' @export
rpoisLocalSC_normal <- nimbleFunction(
  run = function( n = double(0, default = 1),
                  lambda0 = double(0),
                  sigma = double(0),
                  s = double(2),
                  trapCoords = double(2),
                  oper = double(1, default = 1.0),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),
                  indicator = double(1)
  ){
    ## Specify return type
    returnType(double(1))
    if(n!=1){print("rpoisLocalSC_normal only allows n = 1; using n = 1")}
    
    ## Identify dimensions
    numIndividuals <- nimDim(s)[1]
    numTraps <- nimDim(trapCoords)[1]
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma)
    
    ## Calculate trap-specific detection rates
    lambda <- nimNumeric( length = numTraps)
    
    ## Loop over individuals 
    for(i in 1:numIndividuals){
      if(indicator[i] > 0){
        ## Retrieve the habitat grid cell ID for this individual's AC 
        sID <- habitatGrid[trunc(s[i,2]/resizeFactor)+1, trunc(s[i,1]/resizeFactor)+1]
        ## Retrieve the local traps surrounding the selected habitat grid cell
        theseLocalTraps <- localTrapsIndices[sID,1:localTrapsNum[sID]]
        ## Loop over local traps only
        for(r in 1:localTrapsNum[sID]){
          ## Calculate distance between AC and local trap
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[i,1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[i,2], 2)
          ## Calculate Poisson rate at this local trap 
          ## This corresponds to the sum of individual lambdas, which follow half-normal distribution
          lambda[theseLocalTraps[r]] <- lambda[theseLocalTraps[r]] + lambda0 * exp(alpha * d2)
        }#r
      }#if
    }#i
    
    ## Calculate the log-probability of the vector of detections
    out <- nimNumeric(numTraps)
    ## Loop over traps 
    for(j in 1:numTraps){
      out[j] <- rpois(n = 1,
                      lambda = lambda[j]*oper[j])
      
    }#j
    
    ## Output
    return(out)
  })

