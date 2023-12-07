dcustom <- nimbleFunction(
  run = function( x = double(0),           ## Number of wolves detected together 
                  size = double(1),        ## Group sizes
                  p0 = double(0),
                  sigma = double(0),
                  alpha = double(0),
                  s = double(2),
                  trapCoords = double(1),
                  localMatrix = double(2),
                  resizeFactor = double(0, default = 1), 
                  habitatGrid = doble(2),                  
                  indicator = double(1),
                  log = integer(0, default = 0)
  ){
    returnType(double(0))
    
    ##-- Calculate group-specific probability of detection :
    ##-- (with local evaluation)
    numGroups <- length(size)
    p <- nimNumeric(length = numGroups)
    sID.trap <- habitatGrid[trunc(trapCoords[2]/resizeFactor) + 1,
                            trunc(trapCoords[1]/resizeFactor) + 1]
    for(g in 1:numGroups){
      if(indicator[g] > 0){
        ##-- Identify habitat cell for this group AC
        sID.group <- habitatGrid[trunc(s[g,2]/resizeFactor) + 1,
                                 trunc(s[g,1]/resizeFactor) + 1]
        
        if(localMatrix[sID.group,sID.trap] > 0){
          ##-- Distance between group AC and trap coords
          d2 <- (s[g,1]-trapCoords[1])^2 + (s[g,2]-trapCoords[2])^2
          ##-- Group-specific probability of visit:
          p[g] <- p0 * exp(-d2 / (2*sigma^2)) 
        }
      }
    }#g
    
    
    ##-- Derive probability of 0 visit:
    sum.p <- sum(p)
    risk <- 1 - exp(-sum.p)
    p <- risk * p/sum.p
    pNULL <- 1 - risk
    logProb <- dmulti(y, prob = p, size = size, log = TRUE)
    
    
    ##-- if n = 0, return probability of 0 visit :
    if(x == 0){
        if(log == 0) 
          return(pNULL)
        else return(log(pNULL))
      }


    ##-- If x > 0; loop over groups and calculate the probability
    ##-- to detect X individuals simultaneously given the group size 
    ##-- and probability of visit :
    for(g in 1:numGroups){
      if(x <= size[g]){
      thisProb <- dtbinom( x,
                           prob = alpha,
                           size = size[g],
                           a = 1,
                           log = 0)
      totalProb <- totalProb + thisProb * p[g]
      }#if
    }#g
    if(log) return(log(thisProb)) else return(thisProb)
  })

