#' @title NIMBLE Function to calculate AC density.
#'
#' @description
#' \code{calculateDensity} is a NIMBLE function to calculate the number of ACs in each habitat grid cell.
#' 
#' @param s \code{matrix} of dimensions n.individuals * 2 containing individual activity center coordinates (scaled to the habitat). 
#' @param habitatGrid \code{matrix} containing IDs of habitat grid cells.
#' @param indicator \code{vector} of length n.individuals denoting if each individual has to be considered in the density calculation or not (e.g. if the individual is dead)
#' @param numWindows \code{numeric} denoting the number of habitat grid cells.

#' @examples
#' d2[i,1:n.detectors,t] <- calculateDistance(sxy[i,1:2,t], detector.xy[1:n.detectors,1:2,t])

calculateDensity <- nimbleFunction(run = function( s = double(2),
                                                   habitatGrid = double(2),
                                                   indicator = double(1),
                                                   numWindows = double(0)){
  # Return type declaration
  returnType(double(1))
  
  n.individuals <- dim(s)[1]
  
  D <- numeric(length = numWindows, value = 0)
  
  for(i in 1:n.individuals){
    if(indicator[i] == 1){
    windowInd <- habitatGrid[trunc(s[i,2]) + 1, trunc(s[i,1]) + 1]
    D[windowInd] <- D[windowInd] + 1
    }
  }#i
  
   return(D)
})