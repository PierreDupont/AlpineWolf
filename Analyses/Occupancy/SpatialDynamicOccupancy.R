
## Spatial dynamic Occupancy Model
library(nimble)

model <- nimbleCode({
  
  ## Ecological submodel
  
  ## Priors for the variance of CAR model
  psi_tau ~ T(dnorm(0, sd = 10), 0, ) # Truncated normal distribution with mean 0 and SD = 10
  
  phi_tau ~ T(dnorm(0, sd = 10), 0, ) # Truncated normal distribution with mean 0 and SD = 10
  
  rho_tau ~ T(dnorm(0, sd = 10), 0, ) # Truncated normal distribution with mean 0 and SD = 10
  
  ## Spatial random effects - CAR model
  psi_rho[1:nsite] ~ dcar_normal(adj=adj[1:sumneigh],
                                 weights=weights[1:sumneigh],
                                 num=num_adj[1:nsite],
                                 tau=psi_tau,
                                 zero_mean = 1)
  
  phi_rho[1:nsite] ~ dcar_normal(adj=adj[1:sumneigh],
                                 weights=weights[1:sumneigh],
                                 num=num_adj[1:nsite],
                                 tau=psi_tau,
                                 zero_mean = 1)
  
  rho_rho[1:nsite] ~ dcar_normal(adj=adj[1:sumneigh],
                                 weights=weights[1:sumneigh],
                                 num=num_adj[1:nsite],
                                 tau=psi_tau,
                                 zero_mean = 1)
  
  ## Spatial fixed effects 
  psi_alpha <- logit(mean.psi)
  mean.psi ~ dunif(0,1)
  psi_betas[1:n.covs] ~ dnorm(0, 0.01)

  phi_alpha <- logit(mean.phi)
  mean.phi ~ dunif(0,1)
  phi_betas[1:n.covs] ~ dnorm(0, 0.01)
  
  rho_alpha <- logit(mean.rho)
  mean.rho ~ dunif(0,1)
  rho_betas[1:n.covs] ~ dnorm(0, 0.01)
  
  
  for (i in 1:nsite){
    ## Initial occupancy
    logit(psi[i]) <- psi_alpha + psi_betas[1:n.covs]%*%covs[i,1:n.covs] + psi_rho[i]
    z[i,1] ~ dbern(psi[i])
    
    ## Colonization/survival processes
    for (k in 2:nyear){
      logit(rho[i,k-1]) <- rho_alpha + rho_betas[1:n.covs]%*%covs[i,1:n.covs] + rho_rho[i]
      
      logit(phi[i,k-1]) <- phi_alpha + phi_betas[1:n.covs]%*%covs[i,1:n.covs] + phi_rho[i]
      
      z[i,k] ~ dbern(z[i,k-1]*phi[i,k-1] + z[i,k-1]*rho[i,k-1])
    }#k
  }#i
  
  
  
  ## Observational model
  p_int   ~ dnorm(0, 0.01)
  
  for (i in 1:nsite){
    for (k in 1:nyear){
      for (j in 1:n.occ[i,k]){ 
        logit(p[i,j,k]) <- p_int  
        muy[i,j,k] <- z[i,k]*p[i,j,k]
        y[i,j,k] ~ dbern(muy[i,j,k])
      }#j
    }#k
  }#i
  
  # Derived parameters
  for (k in 1:nyear){
    occ[k] <- sum(z[1:nsite,k]) # Occupancy in each year
  }
})
