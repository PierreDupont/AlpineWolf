### load packages and prep workspace ####
library(parallel)
library(coda)
library(nimble)
source("C:/My_documents/AlpineWolf/Analyses/SpatialCount/SC_model.R")
load("C:/My_documents/AlpineWolf/Analyses/SpatialCount/SC_model.R")


nadapt<-1000
niter_upd<-8000
str(datasets_list)
    


### static things to define ####
xlims <- datasets_list[[1]]$simDat$xlim
ylims <- datasets_list[[1]]$simDat$ylim
M <- nrow(datasets_list[[1]]$simSC)
SC_data <- datasets_list[[1]]$simSC
oper <- matrix(1,
                nrow(datasets_list[[1]]$simSC),
                ncol(datasets_list[[1]]$simSC))
inits <- list(sigma=rnorm(1,5),lam0=runif(1), z=rep(1,M),psi=0.5)
nimbledata <- list( n = SC_data)
nimbleconstants <- list(M = M,
                        oper = oper,
                        X = datasets_list[[1]]$traplocs,
                        J = nrow(SC_data),
                        K = ncol(SC_data),
                        xlim = xlims,
                        ylim = ylims, 
                        area=((xlims[2]-xlims[1])*(ylims[2]-ylims[1])))
params <- c("N", "D", "lam0", "sigma","psi","z","s")


model <- nimbleModel( code = caribou_nimble,
                      name = "caribou_nimble",
                      constants = nimbleconstants,
                      data = nimbledata,
                      inits = inits)
Cmodel <- compileNimble(model)
modelConf <- configureMCMC(model)
modelConf$addMonitors(params)
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble(modelMCMC, project = model)
out1 <- runMCMC(CmodelMCMC, niter = nadapt)
#   return(as.mcmc(out1))
# })
nimbleEnd <- Sys.time()
nimbleTime<-nimbleEnd-nimbleStart

### check on the first batch of results, 200% likely need to keep running ####
out.mcmc <- as.mcmc(out1)
traceplot(out.mcmc[, "N"])



##-- Try second version of the SC model
oper2 <- rep(ncol(datasets_list[[1]]$simSC),nrow(datasets_list[[1]]$simSC))

nimbledata2 <- list( n = rowSums(SC_data))

nimbleconstants2 <- list(M = M,
                        oper = oper2,
                        X = datasets_list[[1]]$traplocs,
                        J = nrow(SC_data),
                        xlim = xlims,
                        ylim = ylims, 
                        area=((xlims[2]-xlims[1])*(ylims[2]-ylims[1])))
model <- nimbleModel( code = caribou_nimble2,
                      name = "caribou_nimble2",
                      constants = nimbleconstants2,
                      data = nimbledata2,
                      inits = inits)
Cmodel <- compileNimble(model)
modelConf <- configureMCMC(model)
modelConf$addMonitors(params)
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble(modelMCMC, project = model)
#out2 <- runMCMC(CmodelMCMC, niter = nadapt)
out2 <- runMCMC(CmodelMCMC, niter = 10000)
out2.mcmc <- as.mcmc(out2)
traceplot(out2.mcmc[, "N"])



### If has not converged, continue sampling ####
nimbleStart_upd2 <- Sys.time()
out3 <- clusterEvalQ(cl, {
  out1 <- runMCMC(CmodelMCMC, niter = niter_upd*2)
  return(as.mcmc(out1))
})
nimbleEnd_upd2<- Sys.time()
nimbleTime_upd2<-nimbleEnd_upd2-nimbleStart_upd2

### check on the results ####
out.mcmc.update2 <- as.mcmc(out3)
out.mcmc.bind <- mcmc.list()
for (i in seq_len(nc)) {
  out.mcmc.bind[[i]] <- mcmc(rbind(out.mcmc[[i]],
                                   out.mcmc.update1[[i]],
                                   out.mcmc.update2[[i]]))
}
traceplot(out.mcmc.bind[, "N"])
save.image(file=paste0("outputNimble_",LowStr(commonName),"_",study,"_",year,"_",timeFrame,"_",dataType,".RData"))

### If has not converged, continue sampling ####
nimbleStart_upd3 <- Sys.time()
out4 <- clusterEvalQ(cl, {
  out1 <- runMCMC(CmodelMCMC, niter = niter_upd*4)
  return(as.mcmc(out1))
})
nimbleEnd_upd3<- Sys.time()
nimbleTime_upd3<-nimbleEnd_upd3-nimbleStart_upd3

### check on the results ####
out.mcmc.update3 <- as.mcmc(out4)
out.mcmc.bind <- mcmc.list()
for (i in seq_len(nc)) {
  out.mcmc.bind[[i]] <- mcmc(rbind(out.mcmc[[i]],
                                   out.mcmc.update1[[i]],
                                   out.mcmc.update2[[i]],
                                   out.mcmc.update3[[i]]))
}
traceplot(out.mcmc.bind[, "N"])
save.image(file=paste0("outputNimble_",LowStr(commonName),"_",study,"_",year,"_",timeFrame,"_",dataType,".RData"))


### when finished ###
summary(out.mcmc.bind)#print(jagsfit_all)
gel<-gelman.diag(out.mcmc.bind, autoburnin=FALSE,multivariate=FALSE)
neff<-effectiveSize(out.mcmc.bind)
nits_total<-paste0(dim(out.mcmc.bind[[1]])[1]," iterations kept per each of ",nc ," chains")
nimbleTime_tot<-sum(get(ls(pattern="^nimbleTime")))
save.image(file=paste0("outputNimble_",LowStr(commonName),"_",study,"_",year,"_",timeFrame,"_",dataType,".RData"))
stopCluster(cl)
