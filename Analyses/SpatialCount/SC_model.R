
caribou_nimble<-nimbleCode( {
  sigma ~ dunif(0,10)#dgamma(24,8) # weakly informative prior - 38-619 km2
  lam0 ~ dunif(0,10)
  psi ~ dunif(0,1)#dbeta(1,1)
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    
    for(j in 1:J) {
      distsq[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
      lam[i, j] <- lam0 * exp(-distsq[i,j] / (2*sigma^2)) * z[i] 
    }
  } # End of 1:M
  for(j in 1:J) {
    for(k in 1:K) {
      bigLambda[j, k] <- sum(lam[1:M, j]) * oper[j, k]
      n[j, k] ~ dpois(bigLambda[j, k])
      
    }
  }
  N <- sum(z[1:M])
  D <- N/area
})


caribou_nimble_group<-nimbleCode( {
  sigma ~ dunif(0,10)#dgamma(24,8) # weakly informative prior - 38-619 km2
  p0 ~ dunif(0,1)
  psi ~ dunif(0,1)#dbeta(1,1)
  
  # loop over groups
  for(g in 1:G) {
    z[g] ~ dbern(psi)
    s[g,1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
      upperCoords = upperHabCoords[1:n.habWindows,1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid2[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
    
    for(j in 1:J) {
      distsq[g,j] <- (s[g,1] - X[j,1])^2 + (s[g,2] - X[j,2])^2
      p[g, j] <- p0 * exp(-distsq[g,j] / (2*sigma^2)) * z[g] 
    }
  } # End of 1:G
  for(j in 1:J) {
    for(k in 1:K) {
      n[j, k] ~ dcustomFunction(p[1:G,j],GS[1:G])
      
    }
  }
  N <- sum(z[1:M])
  D <- N/area
})