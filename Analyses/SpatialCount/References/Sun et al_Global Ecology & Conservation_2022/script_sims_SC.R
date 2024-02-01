###fitting SC model
library(parallel)
library(igraph)
library(coda)
library(nimble)

load(file = "simDatasets_Poiss_forSCandSPIM.rda")

prior<-"gamma"
source("CA_nimble_weaklyInformativePriors.R")

M<-400
#xlim
#ylim

buffer<-datasets_list[[1]]$buffer
traplocs<-datasets_list[[1]]$traplocs
X<-as.matrix(traplocs[,-1])
xlim <- c(min(X[, 1]), max(X[, 1])) + c(-buffer, buffer)
ylim <- c(min(X[, 2]), max(X[, 2])) + c(-buffer, buffer)


### static things to define ####

nimbleInits <- list(sigma=5,lam0=0.5, z=rep(1,M),psi=0.5,
                    s=cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2])))
simSC<-datasets_list

nimbleconstants <- list(M = M,
                        oper=matrix(1,nrow=dim(simSC[[1]]$simSC)[1],ncol=dim(simSC[[1]]$simSC)[2]),
                        X=X,J = nrow(simSC[[1]]$simSC),
                        K=ncol(simSC[[1]]$simSC),
                        xlim=xlim, ylim=ylim, 
                        area=((xlim[2]-xlim[1])*(ylim[2]-ylim[1])))
params <- c("N", "D", "lam0", "sigma")#,"psi","z","s")



niter<-50000

itsToDo<-c(1:50)

### set up cluster ####
set.seed(21) #24
nc <- length(itsToDo) # number of sims
cl<-makeCluster(nc,timeout=5184000)

### tell the cluster what it will need, things that dont change ####
clusterExport(cl, c("niter","caribou_nimble","simSC",
                    "nimbleInits", "nimbleconstants", "params", "itsToDo"))

#i think this for loop could be removed if nimbledata was moved to the call above
# and the jth dataset extracted in the clusterEvalQ below, 
#but i inherited this structure when the inits were a function dependent on the set.seed.
for (j in seq_along(cl)) {
  
  #j and k will be the same if itsToDo is sequential from 1
  k<-itsToDo[j]
  nimbledata <- list( n = simSC[[k]]$simSC)
  
  clusterExport(cl[j], c("j","k","nimbledata"))
}

### tell the cluster what to do ####
#out <- 
clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  #testing if this call structure works
  #test<-niter+j
  word<-paste("this was iteration ", k)
  save.image(file=paste0("outputNimble_simSC_",k,".RData"))
  
  nimbleStart <- Sys.time()
  #nimbleOptions(MCMCprogressBar=TRUE,showCompilerOutput=TRUE)
  
  model <- nimbleModel(code = caribou_nimble, name = "caribou_nimble",
                       constants = nimbleconstants, data = nimbledata,
                       inits = nimbleInits)
  Cmodel <- compileNimble(model)
  
  nimbleStart_configure <- Sys.time()
  modelConf <- configureMCMC(model)
  modelConf$addMonitors(params)
  modelMCMC <- buildMCMC(modelConf)
  nimbleEnd_modelBuilding <- Sys.time()
  
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  
  nimbleStart_run <- Sys.time()
  out1 <- runMCMC(CmodelMCMC, niter = niter)
  save.image(file=paste0("outputNimble_simSC_",k,".RData"))
  
  nimbleEnd <- Sys.time()
  
  out.mcmc <- as.mcmc(out1)
  neff<-effectiveSize(out.mcmc)
  
  nimbleTime_total<-nimbleEnd-nimbleStart
  nimbleTime_modelBuilding<-nimbleEnd_modelBuilding-nimbleStart_configure
  nimbleTime_run<-nimbleEnd-nimbleStart_run
  
  save.image(file=paste0("outputNimble_simSC_",k,".RData"))
  
  #return(as.mcmc(out1))
})


stopCluster(cl)

