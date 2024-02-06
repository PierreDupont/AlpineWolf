calcLambda <- nimbleFunction(
  run = function( lambda0 = double(0),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  indicator = double(0), 
                  log = integer(0, default = 0)
  ){
    ## Specify return type
    returnType(double(1))
    
    ## Identify dimensions
    numTraps <- nimDim(trapCoords)[1]
    
    ## Calculate trap-specific detection rates
    lambda <- nimNumeric(length = numTraps)
    
    ## Shortcut if individual is not available for detection
    if (indicator < 1){return(lambda)}
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma) 
    
    ## Loop over traps
    for(j in 1:numTraps){
      d2 <- pow(trapCoords[j,1] - s[1], 2) + pow(trapCoords[j,2] - s[2], 2)
      lambda[j] <- lambda0 * exp(alpha * d2)
    }
    return(lambda)
  })


calcLocalLambda <- nimbleFunction(
  run = function( lambda0 = double(0),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),
                  indicator = double(0), 
                  log = integer(0, default = 0)
  ){
    ## Specify return type
    returnType(double(1))
    
    ## Identify dimensions
    numTraps <- nimDim(trapCoords)[1]
    
    ## Calculate trap-specific detection rates
    lambda <- nimNumeric(length = numTraps)
    
    ## Shortcut if individual is not available for detection
    if (indicator < 1){return(lambda)}
    
    ## Initialize objects
    alpha <- -1.0 / (2.0 * sigma * sigma)
    
    ## Retrieve the habitat grid cell ID for this individual's AC 
    sID <- habitatGrid[trunc(s[i,2]/resizeFactor)+1, trunc(s[i,1]/resizeFactor)+1]
    
    ## Retrieve the local traps surrounding the selected habitat grid cell
    theseLocalTraps <- localTrapsIndices[sID,1:localTrapsNum[sID]]
    
    ## Loop over local traps only
    for(r in 1:localTrapsNum[sID]){
      ## Calculate distance between AC and local trap
      d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
      ## Calculate Poisson rate at this local trap 
      ## This corresponds to the sum of individual lambdas, which follow half-normal distribution
      lambda[theseLocalTraps[r]] <- lambda0 * exp(alpha * d2)
    }#r
    
    return(lambda)
  })