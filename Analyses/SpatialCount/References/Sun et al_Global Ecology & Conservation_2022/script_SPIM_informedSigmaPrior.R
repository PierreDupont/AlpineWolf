###fitting spim model
#packages
library(devtools)
install_github("benaug/SPIM")
library(SPIM)
library(coda)
library(tidyr)
library(ggplot2)


#data and prep
study<-"Algar"
species1<-"caribou"
timeFrame<-"3mos"
dataType<-"all"#grpSize1" #all
year<- 2019
chain<-1

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

load(paste0("data_",species1,"_",study,"_",year,"_3mos_forSPIM.RData"))

SPIM_data<-get(paste0(commonName,"_SPIM_",year,"_",dataType))
SC_camOper<-get(paste0(LowStr(commonName),"_camOp_",year))
table(rowSums(SC_camOper))
dim(SC_camOper)

#traplocs J
head(traplocs)
dim(traplocs)
camOp_focal<-get(paste0(species1,"_camOp_",year))
rownames(camOp_focal) #sites that operated at least one day and thus are included in the traplocs

#trap buffer for state space
buffer<-2.32*3 # 

#obstype  (poisson, bernoulli)
obstype<-"poisson"


#occasions K
ncol(camOp_focal) #number of days in season
colnames(camOp_focal) #dates in season

#Data, keeping only obs that are on the dates and at the sites in camop
y.obs_long<-get(paste0(CapStr(species1),"_SPIM_",year,"_",dataType))
obs_withSitesToKeep<-which(y.obs_long$Deployment.Location.ID%in%rownames(camOp_focal))
obs_withDatesToKeep<-which(as.character(as.Date(y.obs_long$Date_Time.Captured))%in%colnames(camOp_focal))
obs_toKeep<-intersect(obs_withSitesToKeep,obs_withDatesToKeep)
y.obs_long<-y.obs_long[obs_toKeep,]
dim(y.obs_long)

y.obs<-array(0,dim=c(nrow(y.obs_long),nrow(traplocs),ncol(camOp_focal)))
for(i in 1:nrow(y.obs_long)){
  whichInd<-i
  whichLoc<-which(rownames(camOp_focal)%in%as.character(y.obs_long$Deployment.Location.ID[i]))
  whichOcc<-which(colnames(camOp_focal)%in%as.character(as.Date(y.obs_long$Date_Time.Captured[i])))
  print(paste0("Ind ",i," at ",whichLoc, " on ",whichOcc))
  y.obs[whichInd,whichLoc,whichOcc]<-1
}

#catvar values
#change missing sex, which are blanks, to unknown, then 3 and then 0; female is 1, male is 2
#change antler counts of 0 to 999 and missing counts to 0
#add a collar column where 2 is yes, 1 is no, and 0 is missing
y.obs_long$Sex<-as.character(y.obs_long$Sex)
y.obs_long$Sex[c(which(y.obs_long$Sex==""),which(is.na(y.obs_long$Sex)))]<-"unknown"
y.obs_long$Sex[which(y.obs_long$Sex%in%c("Male","M","male"))]<-"Male"
y.obs_long$Sex[which(y.obs_long$Sex%in%c("Female","f","female"))]<-"Female"
y.obs_long$Sex<-as.numeric(as.factor(y.obs_long$Sex))
y.obs_long$Sex[which(y.obs_long$Sex==3)]<-0

y.obs_long$Antler.pts[which(y.obs_long$Antler.pts==0)]<-999
y.obs_long$Antler.pts[which(is.na(y.obs_long$Antler.pts))]<-0

G.obs<-cbind(y.obs_long$Sex,y.obs_long$Antler.pts) # i . # cat vars

#if also wanting to include collar information
y.obs_long$Collar<-as.numeric(y.obs_long$Comments=="has collar")+1
G.obs<-cbind(y.obs_long$Sex,y.obs_long$Antler.pts,y.obs_long$Collar) 
head(G.obs)

#list of catvars
ncat<-ncol(G.obs)
IDcovs<-vector("list",ncat)#Store unique genotypes
for(c in 1:length(IDcovs)){
  if(max(G.obs[,c])==999){ #if working on antler cov and there are missing values, ie 999
    #because otherwise it will create categories all the way up to 999 which you dont want
    IDcovs[[c]]<-c(1:(sort(unique(G.obs[,c]), TRUE)[2]),999) 
  }else{
    IDcovs[[c]]<-1:max(G.obs[,c],na.rm = TRUE)  
  }
  
}
IDcovs
IDlist<-list(ncat=ncat,IDcovs=IDcovs)

#gamma, the initial values for pop level proportions for each catvar
gamma_props<-vector("list",ncat)#Store unique genotypes
for(c in 1:2){
  gamma_props[[c]]<-rep(1/length(IDcovs[[c]]),length(IDcovs[[c]]))
}
if(ncat==3){
  gamma_props[[3]]<-c((100-length(which(G.obs[,3]==2)))/100,length(which(G.obs[,3]==2))/100)
}
gamma_props

#initial value for sigma
sigma<-5

# wrap data int a list
dataSPIM<-list(y.obs=y.obs,G.obs=G.obs,IDlist=IDlist,
               X=traplocs,K=ncol(camOp_focal),
               buff=buffer, obstype=obstype)


#MCMC stuff
niter=500000 #how long to run the MCMC chain. & minutes for 10k to get  
nburn=0 #how much burnin to discard. I always do this afterwards so I can assess convergence better
nthin=1 #do we thin the chain. Only necessary if you need to reduce file size
nswap=nrow(y.obs)/2 #Number of latent IDs to swap on each iteration. Updating half seems to work fine.
#Theoretically a tradeoff between mixing and run time, but mixing is fine with not many swaps.
M=400 #Data augmentation 
inits=list(psi=0.5,lam0=0.5,sigma=sigma,gamma=gamma_props) #initial values. Using simulated values
#initial values. Using simulated values
priors=list(sigma=c(24,8))
proppars=list(lam0=0.05,sigma=0.075,sx=1,sy=1) #MCMC proposal parameters. Tuning these is the most difficult part.
proppars=list(lam0=0.0025,sigma=0.33,sx=1,sy=1)
#shoot for rejection rates between 0.2 and 0.4. 
#If a parameter is not moving, the proposal parameter is too large.
IDup="MH" #Gibbs or metropolis-hastings latent ID update. Must use MH for binomial model.
#Both about equally as fast for Poisson
keepACs=TRUE #Do we store the activity center posteriors and other latent structure including the ID vector?
keepGamma=TRUE #Do we store the gamma posterior?

if(chain==1){
  set.seed(123)
}
if(chain==3){
  set.seed(789)
}
if(chain==2){
  set.seed(456)
}
a=Sys.time()
out<-mcmc.CatSPIM(dataSPIM,tf=camOp_focal,
                  niter=niter,nburn=nburn,nthin=nthin, nswap=nswap,
                  M = M, inits=inits,proppars=proppars,obstype=obstype,priors=priors,
                  IDup=IDup,keepACs=keepACs,keepGamma=keepGamma)
b=Sys.time()
b-a
save.image(file=paste0("outputSPIM_",LowStr(commonName),"_",study,"_",year,"_",timeFrame,"_",dataType,"_chain",chain,".RData"))


###FUN THINGS

#This will get you  acceptance probabilies. Can't change N, n, or psi.
acceptRate<-1-rejectionRate(mcmc(out$out))
acceptRate
#effective sample size;should shoot for at least 400 for N. 
effectiveSize(out$out)


#Let's see what happened!
plot(mcmc(out$out))
plot(mcmc(out$out[1000:niter,]))

#If you kept the gamma posteriors, you can plot those, too.
plot(mcmc(out$gammaOut[[1]]))
plot(mcmc(out$gammaOut[[2]]))
gamma #True values

#This will get you the acceptance probabilities for the x dimension of the activity centers.
#Hard to tune because it will be lower when an activity center has samples allocated
#to it and higher when it does not. I think just make sure it's not under 0.1ish
1-rejectionRate(mcmc(out$sxout))

#Get posterior probability sample x and sample y came from same individual
niters=niter-nburn

check=1 #Which sample to check
storematch=matrix(0,nrow=ncol(out$IDout),ncol=niters)
for(i in 1:niters){
  storematch[,i]=out$IDout[i,check]==out$IDout[i,]
}
rowSums(storematch)/niters
#True IDs stored here
#data$IDtrue


PID<-matrix(0,nrow=ncol(out$IDout),ncol=1+ncol(out$IDout))
PID[,1]<-seq(1,ncol(out$IDout),1)
PID<-as.data.frame(PID)
colnames(PID)<-c("ID",seq(1,ncol(out$IDout),1))
rownames(PID)<-colnames(PID)[-1]
for(z in 1:ncol(out$IDout)){
  print(z)
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

ggplot(PID, aes(x=ID, y=ComparedTo_Ind)) +
  geom_tile(aes(fill = Probability)) +
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title = "posterior probability sample x and sample y came from same individual",
       x="ID",y="ID")

