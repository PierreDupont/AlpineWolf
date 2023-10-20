### load pacakges and prep workspace ####
library(parallel)
library(coda)
library(nimble)
source("SC_model.R")

### dynamic things to define ####
study<-"Rich" #"Algar"
species<-"caribou"
year<- 2018
timeFrame<- "3mos"
dataType<-"sum"

nadapt<-1000
niter_upd<-8000

CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}
LowStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(tolower(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

load(paste0("data_",species,"_",study,"_",year,"_",timeFrame,"_forSC.RData"))
SC_data<-get(paste0(LowStr(commonName),"_SC_",CapStr(dataType),"_",year))
SC_camOper<-get(paste0(LowStr(commonName),"_camOp_",year))
             


### static things to define ####
             
inits <- function()list(sigma=rnorm(1,5),lam0=runif(1), z=rep(1,M),psi=0.5)
nimbledata <- list( n = SC_data)
nimbleconstants <- list(M = M,
                        oper=SC_camOper,
                        X=traplocs,J = nrow(SC_data),
                        K=ncol(SC_data),
                        xlim=xlims, ylim=ylims, 
                        area=((xlims[2]-xlims[1])*(ylims[2]-ylims[1])))
params <- c("N", "D", "lam0", "sigma","psi","z","s")



### set up cluster ####
nimbleStart <- Sys.time()
set.seed(22)
nc <- 3 # number of chains
cl<-makeCluster(nc,timeout=5184000)

### tell the cluster what it will need ####
clusterExport(cl, c("nadapt","niter_upd","caribou_nimble", "inits", "nimbledata", "nimbleconstants", "params"))
for (j in seq_along(cl)) {
  set.seed(j)
  init <- inits()
  clusterExport(cl[j], "init")
}

### tell the cluster what to do ####
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  model <- nimbleModel(code = caribou_nimble, name = "caribou_nimble",
                       constants = nimbleconstants, data = nimbledata,
                       inits = init)
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  modelConf$addMonitors(params)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  out1 <- runMCMC(CmodelMCMC, niter = nadapt)
  return(as.mcmc(out1))
})
nimbleEnd <- Sys.time()
nimbleTime<-nimbleEnd-nimbleStart

### check on the first batch of results, 200% likely need to keep running ####
out.mcmc <- as.mcmc(out)
traceplot(out.mcmc[, "N"])

### If has not converged, continue sampling ####
nimbleStart_upd1 <- Sys.time()
out2 <- clusterEvalQ(cl, {
  out1 <- runMCMC(CmodelMCMC, niter = niter_upd)
  return(as.mcmc(out1))
})
nimbleEnd_upd1<- Sys.time()
nimbleTime_upd1<-nimbleEnd_upd1-nimbleStart_upd1

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
