#1.	Generate population and SCR encounter histories 
#2.	Simulate encounter histories with sampling design
#2a.	Reduce to SC data by removing identities and summing 
#2b.	Reduce to SPIM data by keeping only some of the ID covs
#3.	Run models 
#3a.   SC model  
#3b.   SPIM models

library(extraDistr)
library(dplyr)
library(ggplot2)


# adapted sim function from oSCR, from Cat Sun
# changed to also output activity centers
sim.SGM <- function(G = 10,                ## Number of groups
                    K = 20,                ## Number of detection occasions
                    p0 = 0.75,             ## Baseline probability of visit 
                    sigma = 0.5,           ## Scale parameter
                    lambda = 3,            ## Mean group size in the population
                    alpha = 0.75,          ## Cohesion of the group (probability to detect X ids together)
                    ssRes = 0.5,           ## Habitat resolution
                    traps_dim = c(5,5),    ## Detector grid size
                    trapCoords = NULL, 
                    buffer = NULL){
  require(extraDistr)
  ##-- Set up trap locations
  if (is.null(trapCoords)) {
    trapCoords <- expand.grid(X = seq(1, traps_dim[1], by = 1), 
                            Y = seq(1, traps_dim[2], by = 1))
  }
  ntraps <- nrow(trapCoords)
  
  ##-- Calculate buffer size
  if (is.null(buffer)) { buffer <- sigma * 3}
  
  ##-- Set up habitat boundaries of the rectangular habitat
  Xl <- min(trapCoords[, 1] - buffer)
  Xu <- max(trapCoords[, 1] + buffer)
  Yl <- min(trapCoords[, 2] - buffer)
  Yu <- max(trapCoords[, 2] + buffer)
  
  ##-- Sample group activity centers
  sx <- runif(G, Xl, Xu)
  sy <- runif(G, Yl, Yu)
  S <- cbind(sx, sy)
  
  ##-- Sample group sizes
  GS <- rtpois( n = G,
                lambda = lambda,
                a = 0)
  
  ##-- Calculate average group density in the population
  Dens <- G/((Xu - Xl) * (Yu - Yl))
  
  ##-- Function to calculate distances
  e2dist<-function (x, y){
    if (!is.matrix(x)) 
      x <- as.matrix(x)
    if (!is.matrix(y)) 
      y <- as.matrix(y)
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
  ##-- Calculate distance between all traps and activity centers
  D <- e2dist(S, trapCoords)
  
  ##-- Calculate alpha for faster calculation
  alpha1 <- 1/(2 * sigma * sigma)
  
  ##-- Calculate group- and trap-specific visit probabilities
  p <- p0 * exp(-alpha1 * D * D)
  
  ##-- Set-up observation matrix y:
  ##-- the dimensions are ntraps*K and it contains the number of wolves 
  ##-- photographed together for each trap and occasion
  Y <- matrix(0, nrow = ntraps, ncol = K)
  dimnames(Y) <- list("traps" = paste("trap", 1:nrow(Y)),
                      "occasions" = 1:ncol(Y))
  ##-- Loop over occasions
  for(k in 1:K){
    
  ##-- Sample group visit (or absence of visit)
  V <- apply(p, 2, function(thisCol){
    sum.p <- sum(thisCol)
    risk <- 1 - exp(-sum.p)
    thisCol <- risk * thisCol/sum.p
    newP <- c(thisCol, 1 - risk) 
    out <- rmulti(n = 1,
                  size = 1,
                  prob = newP)
    return(out)
  })
  

  
  ##-- Sample number of wolves detected together for each visit
  for(j in 1:ntraps){
    ##-- Identify which group visited detector j (if any)
    whichGroup <- which(V[1:G,j] > 0)
    ##-- Sample number of wolves detected together
    if(length(whichGroup)>0){
      Y[j,k] <- rtbinom(n = 1,
                        size = GS[whichGroup],
                        prob = alpha,
                        a = 0)
    }#if
  }#j
}#k
  
  ##-- Output
  return(list( "Y" = Y,
               "GS" = GS, 
               "S" = S,
               "trapCoords" = trapCoords,
               "habitatBoundaries" = list("x" = c(Xl,Xu),
                                          "y" = c(Yl,Yu)),
               "visits" = V))
}



###1 - generate population and simulate detections ####
testSim <- sim.SGM()
