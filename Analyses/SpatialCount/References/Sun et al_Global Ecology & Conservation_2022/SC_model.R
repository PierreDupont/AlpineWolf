
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
