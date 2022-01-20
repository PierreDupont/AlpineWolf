pvaluefn <- function(repData = NULL,# list (niter) # components: matrix (M x ntrap) 
                     obsData = NULL,# M x ntrap
                     reploglik = NULL, # log-likelihood # list (niter) # components: vectors (M x 1)
                     obsloglik = NULL, # log-likelihood # list (niter) # components: vectors (M x 1)
                     detprob = NULL, # list (niter) # components: matrix (M x ntrap) 
                     zchain = NULL, # list (niter) # components: vectors (M x 1) 
                     ntrial = 1, # no. of occasions or trial
                     metric = 'FT') # distance_metric = c('Freeman_Tukey', 'omnibus', 'LR', 'LR2', 'deviance', 'sqrt', 'Freeman-Tukey2')
  
{
  statnames <- c('FT', 'FT_detector', 'FT_individual', 'Omnibus', 
                 'FT_z', 'FT_detector_z', 'FT_individual_z', 'Omnibus_z', 
                 'LR', 'Deviance', 'sqrt', 'sqrt_z')
  if(sum(metric %in% 'ALL') == 0) metric <- intersect(statnames, metric)
  
  Tstat.sim <- Tstat.obs <- NULL
  pvalue.out <- c()
  if(sum(metric %in% c('FT', 'ALL')) > 0){
    disfn <- function(y1, pr){ sum((sqrt(y1) - sqrt(pr*ntrial))^2) }
    Tstat.sim.FT <- mapply(disfn, y1 = repData, pr = detprob, SIMPLIFY = T) # niter x 1
    Tstat.obs.FT <- unlist(lapply(detprob, disfn, y1 = obsData)) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.FT)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.FT)
    if(sum(metric %in% c('FT')) > 0){
      pvalue.FT <- mean(Tstat.sim.FT>Tstat.obs.FT)  
      pvalue.out <- c(pvalue.out, pvalue.FT)
    }
  }
  if(sum(metric %in% c('FT_detector', 'ALL')) > 0){
    disfn <- function(y1, pr){ sum((sqrt(colSums(y1)) - sqrt(colSums(pr)*ntrial))^2) }
    Tstat.sim.FT.det <- mapply(disfn, y1 = repData, pr = detprob, SIMPLIFY = T) # niter x 1
    Tstat.obs.FT.det <- unlist(lapply(detprob, disfn, y1 = obsData)) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.FT.det)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.FT.det)
    if(sum(metric %in% c('FT_detector')) > 0){
      pvalue.FT.det <- mean(Tstat.sim.FT.det>Tstat.obs.FT.det)  
      pvalue.out <- c(pvalue.out, pvalue.FT.det)
    }
  }
  if(sum(metric %in% c('FT_individual', 'ALL')) > 0){
    disfn = function(y1, pr){ sum((sqrt(rowSums(y1)) - sqrt(rowSums(pr)*ntrial))^2) }
    Tstat.sim.FT.ind <- mapply(disfn, y1 = repData, pr = detprob, SIMPLIFY = T) # niter x 1
    Tstat.obs.FT.ind <- unlist(lapply(detprob, disfn, y1 = obsData)) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.FT.ind)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.FT.ind)
    if(sum(metric %in% c('FT_individual')) > 0){
      pvalue.FT.ind <- mean(Tstat.sim.FT.ind>Tstat.obs.FT.ind)  
      pvalue.out <- c(pvalue.out, pvalue.FT.ind)
    }
  }
  if(sum(metric %in% c('Omnibus', 'ALL')) > 0){
    disfn = function(y1, pr){ sum(((y1 - pr*ntrial)^2)/(pr*ntrial)) }
    Tstat.sim.omni <- mapply(disfn, y1 = repData, pr = detprob, SIMPLIFY = T) # niter x 1
    Tstat.obs.omni <- unlist(lapply(detprob, disfn, y1 = obsData)) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.omni)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.omni)
    if(sum(metric %in% c('Omnibus')) > 0){
      pvalue.omni <- mean(Tstat.sim.omni>Tstat.obs.omni)  
      pvalue.out <- c(pvalue.out, pvalue.omni)
    }
  }
  if(sum(metric %in% c('FT_z', 'ALL')) > 0){
    disfn <- function(y1, pr, zz){ sum(rowSums((sqrt(y1) - sqrt(pr*ntrial))^2)*zz) }
    Tstat.sim.FT.z <- mapply(disfn, y1 = repData, pr = detprob, zz = zchain, SIMPLIFY = T) # niter x 1
    disfn2 <- function(pr, zz){ sum(rowSums((sqrt(obsData) - sqrt(pr*ntrial))^2)*zz) }
    Tstat.obs.FT.z <- mapply(disfn2, pr = detprob, zz = zchain, SIMPLIFY = T) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.FT.z)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.FT.z)
    if(sum(metric %in% c('FT_z')) > 0){
      pvalue.FT.z <- mean(Tstat.sim.FT.z>Tstat.obs.FT.z)  
      pvalue.out <- c(pvalue.out, pvalue.FT.z)
    }
  }
  if(sum(metric %in% c('FT_detector_z', 'ALL')) > 0){
    disfn <- function(y1, pr, zz){ sum((sqrt(colSums(y1*zz)) - sqrt(colSums(pr*zz)*ntrial))^2) }
    Tstat.sim.FT.det.z <- mapply(disfn, y1 = repData, pr = detprob, zz = zchain, SIMPLIFY = T) # niter x 1
    disfn2 <- function(pr, zz){ sum((sqrt(colSums(obsData*zz)) - sqrt(colSums(pr*zz)*ntrial))^2) }
    Tstat.obs.FT.det.z <- mapply(disfn2, pr = detprob, zz = zchain, SIMPLIFY = T) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.FT.det.z)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.FT.det.z)
    if(sum(metric %in% c('FT_detector_z')) > 0){
      pvalue.FT.det.z <- mean(Tstat.sim.FT.det.z>Tstat.obs.FT.det.z)  
      pvalue.out <- c(pvalue.out, pvalue.FT.det.z)
    }
  }
  if(sum(metric %in% c('FT_individual_z', 'ALL')) > 0){
    disfn = function(y1, pr, zz){ sum(((sqrt(rowSums(y1)) - sqrt(rowSums(pr)*ntrial))^2)*zz) }
    Tstat.sim.FT.ind.z <- mapply(disfn, y1 = repData, pr = detprob, zz = zchain, SIMPLIFY = T) # niter x 1
    disfn2 = function(pr, zz){ sum(((sqrt(rowSums(obsData)) - sqrt(rowSums(pr)*ntrial))^2)*zz) }
    Tstat.obs.FT.ind.z <- mapply(disfn2, pr = detprob, zz = zchain, SIMPLIFY = T) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.FT.ind.z)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.FT.ind.z)
    if(sum(metric %in% c('FT_individual_z')) > 0){
      pvalue.FT.ind.z <- mean(Tstat.sim.FT.ind.z>Tstat.obs.FT.ind.z)  
      pvalue.out <- c(pvalue.out, pvalue.FT.ind.z)
    }
  }
  if(sum(metric %in% c('Omnibus_z', 'ALL')) > 0){
    disfn = function(y1, pr, zz){ sum(rowSums(((y1 - pr*ntrial)^2)/(pr*ntrial))*zz) }
    Tstat.sim.omni.z <- mapply(disfn, y1 = repData, pr = detprob, zz = zchain, SIMPLIFY = T) # niter x 1
    disfn2 = function(pr, zz){ sum(rowSums(((obsData - pr*ntrial)^2)/(pr*ntrial))*zz) }
    Tstat.obs.omni.z <- mapply(disfn2, pr = detprob, zz = zchain, SIMPLIFY = T) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.omni.z)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.omni.z)
    if(sum(metric %in% c('Omnibus_z')) > 0){
      pvalue.omni.z <- mean(Tstat.sim.omni.z>Tstat.obs.omni.z)  
      pvalue.out <- c(pvalue.out, pvalue.omni.z)
    }
  }
  
  if(sum(metric %in% c('LR', 'ALL')) > 0){
    disfn = function(y1, pr){ 2*sum(y1*(log(y1+0.0001) - log(pr*ntrial))) }
    Tstat.sim.LR <- mapply(disfn, y1 = repData, pr = detprob, SIMPLIFY = T) # niter x 1
    Tstat.obs.LR <- unlist(lapply(detprob, disfn, y1 = obsData)) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.LR)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.LR)
    if(sum(metric %in% c('LR')) > 0){
      pvalue.LR <- mean(Tstat.sim.LR>Tstat.obs.LR)  
      pvalue.out <- c(pvalue.out, pvalue.LR)
    }
  }
  if(sum(metric %in% c('Deviance', 'ALL')) > 0){
    Tstat.sim.dev <- -2*unlist(lapply(reploglik, sum))
    Tstat.obs.dev <- -2*unlist(lapply(obsloglik, sum))
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.dev)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.dev)
    if(sum(metric %in% c('Deviance')) > 0){
      pvalue.dev <- mean(Tstat.sim.dev>Tstat.obs.dev)  
      pvalue.out <- c(pvalue.out, pvalue.dev)
    }
  }
  if(sum(metric %in% c('sqrt', 'ALL')) > 0){
    disfn = function(y1){ sum(sqrt(rowSums(y1))) }
    Tstat.sim.sqrt <- lapply(repData, disfn)
    Tstat.obs.sqrt <- disfn(obsData)
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.sqrt)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.sqrt)
    if(sum(metric %in% c('sqrt')) > 0){
      pvalue.sqrt <- mean(Tstat.sim.sqrt>Tstat.obs.sqrt)  
      pvalue.out <- c(pvalue.out, pvalue.sqrt)
    }
  }
  if(sum(metric %in% c('sqrt_z', 'ALL')) > 0){
    disfn = function(y1,zz){ sum(sqrt(rowSums(y1)*zz)) }
    disfn2 = function(zz){ sum(sqrt(rowSums(obsData)*zz)) }
    Tstat.sim.sqrt.z <- mapply(disfn, y1 = repData, zz = zchain, SIMPLIFY = T) # niter x 1
    Tstat.obs.sqrt.z <- simplify2array(lapply(zchain, disfn2)) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.sqrt.z)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.sqrt.z)
    if(sum(metric %in% c('sqrt_z')) > 0){
      pvalue.sqrt.z <- mean(Tstat.sim.sqrt.z>Tstat.obs.sqrt.z)  
      pvalue.out <- c(pvalue.out, pvalue.sqrt.z)
    }
  }
  #=========================================================
  if(sum(metric %in% 'ALL') == 0 && sum(metric %in% statnames) < length(statnames)) { names(pvalue.out) <- metric }
  if(metric == 'ALL' || sum(metric %in% statnames) == length(statnames))
  {
    pvalue.out <- mapply(function(tsim,tobs){mean(tsim>tobs)}, tsim = split2list(Tstat.sim, 2), tobs = split2list(Tstat.obs, 2), SIMPLIFY = T)
    names(pvalue.out) <- statnames
  }
  out <- list(Tstat.sim = Tstat.sim, Tstat.obs = Tstat.obs, pvalue = pvalue.out)
}# pvaluefn function
