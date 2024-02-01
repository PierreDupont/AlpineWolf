# dcustom(
#   x = model$n[1,1]
#          ,
#          size = model$groupSize
#          ,
#          p0 = model$p0
#          ,
#          sigma = model$sigma
#          ,
#          alpha = model$alpha
#          ,
#          s = model$s
#          ,
#          trapCoords = nimConstants$trapCoords[1, ]
#          ,
#          indicator = model$z
#          ,
#          numGroups = nimConstants$G
#          ,
#          log = 0
#          )



dcustom <- nimbleFunction(
  run = function( x = double(0),           ## Number of wolves detected together 
                  size = double(1),        ## Group sizes
                  p0 = double(0),
                  sigma = double(0),
                  alpha = double(0),
                  s = double(2),
                  trapCoords = double(1),
                  indicator = double(1),
                  numGroups = double(0),
                  log = integer(0, default = 0)
  ){
    returnType(double(0))
    
    alpha2 <- -1/(2*sigma^2)
    
    # ##-- Calculate group-specific probability of detection :'  
    # p <- rep(0, numGroups)
    # for(g in 1:numGroups){
    #   if(indicator[g] > 0){
    #     ##-- Distance between group AC and trap coords
    #     d2 <- pow(s[g,1]-trapCoords[1],2) + pow(s[g,2]-trapCoords[2],2)
    #     ##-- Group-specific probability of visit:
    #     p[g] <- p0 * exp(d2 * alpha2) 
    #   }
    # }#g
    # 
    ##-- Calculate group-specific probability of detection :
    p <- rep(0, numGroups)
    for(g in 1:numGroups){
        ##-- Distance between group AC and trap coords
        d2 <- pow(s[g,1]-trapCoords[1],2) + pow(s[g,2]-trapCoords[2],2)
        ##-- Group-specific probability of visit:
        p[g] <- as.numeric(p0 * exp(d2 * alpha2) * indicator[g])
    }#g
    
    
    ##-- Derive probability of 0 visit:
    sum.p <- sum(p[1:numGroups])
    risk <- 1 - exp(-sum.p)
    p <- risk * p/sum.p
    pNULL <- 1 - risk
    #logProb <- dmulti(y, prob = p, size = size, log = TRUE)
    
    
    ##-- if n = 0, return probability of 0 visit :
    if(x == 0){
        if(log == 0) {
          return(pNULL) 
        } else {
          return(log(pNULL))
        }
      }


    ##-- If x > 0; loop over groups and calculate the probability
    ##-- to detect X individuals simultaneously given the group size 
    ##-- and probability of visit :
    totalProb <- 0
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





dcustomLocal <- nimbleFunction(
  run = function( x = double(0),           ## Number of wolves detected together 
                  size = double(1),        ## Group sizes
                  p0 = double(0),
                  sigma = double(0),
                  alpha = double(0),
                  s = double(2),
                  trapCoords = double(1),
                  localMatrix = double(2),
                  resizeFactor = double(0, default = 1), 
                  habitatGrid = double(2),                  
                  indicator = double(1),
                  numGroups = double(0),
                  log = integer(0, default = 0)
  ){
    returnType(double(0))
    
    ##-- Calculate group-specific probability of detection :
    ##-- (with local evaluation)
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
    #logProb <- dmulti(y, prob = p, size = size, log = TRUE)
    
    
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





# rcustom <- nimbleFunction(
#   run = function( n = double(0),           ## Number of wolves detected together 
#                   size = double(1),        ## Group sizes
#                   p0 = double(0),
#                   sigma = double(0),
#                   alpha = double(0),
#                   s = double(2),
#                   trapCoords = double(1),
#                   indicator = double(1),
#                   log = integer(0, default = 0)
#   ){
#     returnType(double(0))
#     
#     ##-- Calculate group-specific probability of detection :
#     numGroups <- length(size)
#     groupDet <- nimNumeric(length = numGroups)
#     for(g in 1:numGroups){
#       if(indicator[g] > 0){
#         ##-- Distance between group AC and trap coords
#         d2 <- (s[g,1]-trapCoords[1])^2 + (s[g,2]-trapCoords[2])^2
#         ##-- Group-specific probability of visit:
#         p <- p0 * exp(-d2 / (2*sigma^2)) 
#         ##-- Sample group detection
#         groupDet[g] <- rbinom( 1,
#                                size = size[g],
#                                prob = p)
#       }
#     }#g
#     
#     # ##-- Derive probability of 0 visit:
#     # #sum.p <- sum(p)
#     # risk <- 1 - exp(-sum.p)
#     # newP <- risk * p/sum.p
#     # pNULL <- 1 - risk
#     # newP <- c(newP, pNULL)
#     #
#     # ##-- Sample group visit
#     # groupDet <- nimNumeric(length = numGroups)
#     # for(g in numGroups){
#     #   groupDet[g] <- rbinom()
#     # }
#     # thisGroup <- rcat(newP)
#     # 
#     # ##-- If thisGroup > numGroups; no visit
#     # ##-- else, sample number of individuals detectedfrom the group considered.
#     # if(thisGroup > numGroups){
#     #   thisDet <- 0    
#     # } else {
#     #   thisDet <- rbinom()   
#     # }
#     
#     ##-- Output
#     thisDet <- max(groupDet)
#     return(thisDet) 
#   })
# 
# rcustomLocal <- nimbleFunction(
#   run = function( x = double(0),           
#                   size = double(1),        
#                   p0 = double(0),
#                   sigma = double(0),
#                   alpha = double(0),
#                   s = double(2),
#                   trapCoords = double(1),
#                   localMatrix = double(2),
#                   resizeFactor = double(0, default = 1), 
#                   habitatGrid = double(2),                  
#                   indicator = double(1),
#                   log = integer(0, default = 0)
#   ){
#     returnType(double(0))
#     
#     ##-- Calculate group-specific probability of detection :
#     ##-- (with local evaluation)
#     numGroups <- length(size)
#     groupDet <- nimNumeric(length = numGroups)
#     sID.trap <- habitatGrid[trunc(trapCoords[2]/resizeFactor) + 1,
#                             trunc(trapCoords[1]/resizeFactor) + 1]
#     for(g in 1:numGroups){
#       if(indicator[g] > 0){
#         ##-- Identify habitat cell for this group AC
#         sID.group <- habitatGrid[trunc(s[g,2]/resizeFactor) + 1,
#                                  trunc(s[g,1]/resizeFactor) + 1]
#         
#         if(localMatrix[sID.group,sID.trap] > 0){
#           ##-- Distance between group AC and trap coords
#           d2 <- (s[g,1]-trapCoords[1])^2 + (s[g,2]-trapCoords[2])^2
#           ##-- Group-specific probability of visit:
#           p <- p0 * exp(-d2 / (2*sigma^2)) 
#           ##-- Sample group detection
#           groupDet[g] <- rbinom( 1,
#                                  size = size[g],
#                                  prob = p)
#         }
#       }
#     }#g
#     
#     ##-- Output
#     thisDet <- max(groupDet)
#     return(thisDet) 
#   })