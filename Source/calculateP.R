## This function calculates the probability of detection of a group using a half-normal detection function
## It uses an indicator to make calculation faster (skips calculations for augmented groups)
calculateP <- nimbleFunction(
  run = function( p0 = double(0),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  indicator = double(0)){
    ## Specify return type
    returnType(double(1))
    
    ## Identify dimensions
    numTraps <- nimDim(trapCoords)[1]
    
    ## Calculate trap-specific detection rates
    p <- nimNumeric(length = numTraps)
    
    ## Shortcut if individual is not available for detection
    if (indicator < 1){return(p)}
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma)
    
    ## Loop over traps
    for(j in 1:numTraps){
      d2 <- pow(trapCoords[j,1] - s[1], 2) + pow(trapCoords[j,2] - s[2], 2)
      p[j] <- p0 * exp(alpha * d2)
    }#j
    return(p)
  })
