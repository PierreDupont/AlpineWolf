###fitting spim model
#packages
library(devtools)
Sys.setenv(TZ='UTC')
library(SPIM)
library(coda)
library(tidyr)
library(ggplot2)

set.seed(123)
setwd("E:/SIMS")

load(file = "simDatasets_Poiss_forSCandSPIM.rda")

whichDataVersion<-"complete" #complete or missing 
whichIDs<-"3ID" #3ID or 2ID . ie. with or without Cllars?
whichdataSet<-which(names(datasets_list[[1]])==paste0("dataSPIM_",whichDataVersion,"_",whichIDs))
whichGammaPriors<-which(names(datasets_list[[1]])==paste0("gamma_props_priors_",whichDataVersion,"_",whichIDs))

names(datasets_list[[1]])

itsForLoop<-c(1:5)
for(itNum in itsForLoop){
  print(itNum)
  whichIt<-itNum
  
  dataSPIM<-datasets_list[[itNum]][[whichdataSet]]
  dataSPIM$obstype<-"poisson"
  y.obs<-dataSPIM$y.obs
  gamma<-datasets_list[[itNum]][whichGammaPriors][[1]]

  if(ncol(dataSPIM$X)==3){
    dataSPIM$X<-dataSPIM$X[,-1]
  }
  
  
  #MCMC stuff
  niter=50000#50000#50000 #how long to run the MCMC chain. 8 minutes for 5k to get  
  nburn=0 #how much burnin to discard. I always do this afterwards so I can assess convergence better
  nthin=1 #do we thin the chain. Only necessary if you need to reduce file size
  nswap=nrow(y.obs)/2 #Number of latent IDs to swap on each iteration. Updating half seems to work fine.
  #Theoretically a tradeoff between mixing and run time, but mixing is fine with not many swaps.
  M=400 #Data augmentation level
  inits=list(psi=0.5,lam0=0.5,sigma=rnorm(1,5))
  inits=list(psi=0.5,lam0=0.5,
             sigma=5,
             gamma=gamma) #initial values. Using simulated values
  #initial values. Using simulated values
  priors=list(sigma=c(24,8))
  
  #MCMC proposal parameters. Tuning these is the most difficult part.
  #shoot for rejection rates between 0.2 and 0.4; acceptance between 0.6-0.8 
  #If a parameter is not moving, the proposal parameter is too large.
  proppars=list(lam0=0.05,sigma=0.075,sx=1,sy=1) 
  proppars=list(lam0=0.002,sigma=0.15,sx=1,sy=1)
  
  IDup="Gibbs" #Gibbs or metropolis-hastings latent ID update. Must use MH for binomial model.
  #Both about equally as fast for Poisson
  
  keepACs=TRUE #Do we store the activity center posteriors and other latent structure including the ID vector?
  keepGamma=TRUE #Do we store the gamma posterior?
  
  
  a=Sys.time()
  out<-mcmc.CatSPIM(dataSPIM,
                    niter=niter,nburn=nburn,nthin=nthin, nswap=nswap,
                    M = M, inits=inits,proppars=proppars,priors=priors,#obstype=obstype,
                    IDup=IDup,keepACs=keepACs,keepGamma=keepGamma)
  b=Sys.time()
  b-a
  
  save.image(file=paste0("outputSimSPIM_",whichDataVersion,"_",whichIDs,"_",whichIt,".RData"))
  
  
  ###FUN THINGS
  
  #This will get you  acceptance probabilies. Can't change N, n, or psi.
  acceptRate<-1-rejectionRate(mcmc(out$out))
  acceptRate
  #effective sample size;should shoot for at least 400 for N. 
  effectiveSize(out$out)
  
  
  #Get posterior probability sample x and sample y came from same individual
  niters=niter-nburn
  
  check=1 #Which sample to check
  storematch=matrix(0,nrow=ncol(out$IDout),ncol=niters)
  for(i in 1:niters){
    storematch[,i]=out$IDout[i,check]==out$IDout[i,]
  }
  rowSums(storematch)/niters
  
  PID<-matrix(0,nrow=ncol(out$IDout),ncol=1+ncol(out$IDout))
  PID[,1]<-seq(1,ncol(out$IDout),1)
  PID<-as.data.frame(PID)
  colnames(PID)<-c("ID",seq(1,ncol(out$IDout),1))
  rownames(PID)<-colnames(PID)[-1]
  for(z in 1:ncol(out$IDout)){
    #print(z)
    check=z #Which sample to check
    storematch=matrix(0,nrow=ncol(out$IDout),ncol=niters)
    for(i in 1:niters){
      storematch[,i]=out$IDout[i,check]==out$IDout[i,]
    }
    PID[z,2:ncol(PID)]<-rowSums(storematch)/niters
  }
  PID<-PID%>%
    pivot_longer(-ID,names_to="ComparedTo_Ind",values_to="Probability")
  PID$ComparedTo_Ind<-as.numeric(PID$ComparedTo_Ind)
  
  save.image(file=paste0("outputSimSPIM_",whichDataVersion,"_",whichIDs,"_",whichIt,".RData"))
  
 
}


