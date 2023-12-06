#1.	Generate population and SCR encounter histories using estimated SPIM parameters
#2.	Simulate encounter histories with sampling design
#2a.	Reduce to SC data by removing identities and summing 
#2b.	Reduce to SPIM data by keeping only some of the ID covs
#3.	Run models 
#3a.   SC model  
#3b.   SPIM models

library(extraDistr)
library(dplyr)
library(ggplot2)

# adapted sim function from oSCR
# changed to also output activity centers
sim.SCR<-function (N = 100,
                   K = 20, 
                   p0 = 0.5,
                   sigma = 0.5, 
                   ssRes = 0.5, 
                   traps_dim = c(5, 5),
                   traplocs = NULL, 
                   buffer = NULL){
  ## Set up trap locations
  if (is.null(traplocs)) {
    traplocs <- expand.grid(X = seq(1, traps_dim[1], by = 1), 
                            Y = seq(1, traps_dim[2], by = 1))
  }
  ntraps <- nrow(traplocs)
  
  ## calculate buffer size
  if (is.null(buffer)) {
    buffer <- sigma * 4
  }
  ## set up habitat boundaries
  Xl <- min(traplocs[, 1] - buffer)
  Xu <- max(traplocs[, 1] + buffer)
  Yl <- min(traplocs[, 2] - buffer)
  Yu <- max(traplocs[, 2] + buffer)
  
  ## Sample group activity centers
  sx <- runif(N, Xl, Xu)
  sy <- runif(N, Yl, Yu)
  S <- cbind(sx, sy)
  
  ## calculate density
  Dens <- N/((Xu - Xl) * (Yu - Yl))
  
  ## Calculate distance between traps and activity centers
  D <- e2dist(S, traplocs)
  
  ## calculate alpha for faster calculation
  alpha1 <- 1/(2 * sigma * sigma)
  
  ## Calculate group- and trap-specific visit probabilities
  probcap <- p0 * exp(-alpha1 * D * D)
  
  ## Derive probability of 0 visit
  pNoGroup <- 
  
  Y <- matrix(NA, nrow = ntraps, ncol = K)
  
  V <- matrix(NA, nrow = ntraps, ncol = K)
  
  for(j in 1:ntraps){
    for(k in 1:K){
      V[j,k] <- rcat()
      
    }
    
  }
  
  for (i in 1:nrow(Y)) {
      Y[i, ] <- rbinom(ntraps, K, probcap[i, ])
  }
  dimnames(Y) <- list(1:nrow(Y), paste("trap", 1:ncol(Y), 
                                       sep = ""))
  if (array3d) {
    Y <- array(NA, dim = c(N, ntraps, K))
    for (i in 1:nrow(Y)) {
      for (j in 1:ntraps) {
        if (encmod == "B") {
          Y[i, j, 1:K] <- rbinom(K, 1, probcap[i, j])
        }
        else if (encmod == "P") {
          Y[i, j, 1:K] <- rpois(K, probcap[i, j])
        }
      }
    }
    if (discard0) {
      Y2d <- apply(Y, c(1, 2), sum)
      ncaps <- apply(Y2d, 1, sum)
      Y <- Y[ncaps > 0, , ]
    }
  }
  ss <- expand.grid(X = seq(Xl + ssRes/2, Xu - ssRes/2, ssRes), 
                    Y = seq(Yl + ssRes/2, Yu - ssRes/2, ssRes))
  list(Y = Y, traplocs = traplocs, xlim = c(Xl, Xu), ylim = c(Yl, Yu),
       N = N, alpha0 = alpha0, alpha1 = alpha1, sigma = sigma, S=S,
       K = K, Dens = Dens, ss = ss, n0 = N - nrow(Y))
}

e2dist<-function (x, y){
  if (!is.matrix(x)) 
    x <- as.matrix(x)
  if (!is.matrix(y)) 
    y <- as.matrix(y)
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

###1 - generate population and simulate detections ####

#N median/mode from 
N<-137 

#sigma median/mode from posterior informed SPIM
sigma<-2.1

#median population proportions of ID covariate values
#from 'outputSPIM_caribou_Algar_2019_all_wCollars_postProcessed_Gammas.RData'
antlers_prop<-c(0.015,0.180,0.036,0.015,0.064,0.039,0.049,0.052,0.090,0.078,0.053,0.058,0.055,0.015,0.015,0.036,0.037,0.004)
sex_prop<-c(0.27,0.73) #  F/M
collar_prop<-c(0.96,0.04) #no collar/ collar

#individuals
pop<-data.frame(Ind=seq(1,N,1),
                Antlers=rcat(N,antlers_prop),
                Sex=rcat(N,sex_prop),
                Collar=rcat(N,collar_prop))
#change the largest antler category to 999, to represent 0 antlers
#because 0 will later be used to represent missing covariate values.
pop$Antlers[which(pop$Antlers==length(antlers_prop))]<-999 

#Algar traplocs, which will also define statespace with buffer
traplocs<-read.csv("Station_Covariates_Feb_2020_STANDARD.csv", header=TRUE)
traplocs<-traplocs[,which(colnames(traplocs)%in%c("Deployment.Location.ID","utmE","utmN"))]
traplocs[,-1]<-traplocs[,-1]/1000 #km
head(traplocs)

buffer<-7 #km

#daily sampling occasions : August, September, October
#assumes perfect operation
K_days<-31+30+31 #days
K_halfhrs<-K_days*24/0.5 #half hours

#detection probability to get the same amount of detections as in observed data
#need to play around with this
alpha0<--5.5
plogis(alpha0)


hist_simSC<-data.frame(Count=c(0,1,2,3,4,5,6,7,8,9,10),
                       n=0)
hist_simSPIM<-c()
hasCollars<-c()

#create datasets in a loop.
datasets_list<-list()
for(it in 1:400){
  simDat<-sim.SCR(N=N, K=K_days,alpha0=alpha0,sigma=sigma,
                  array3d = TRUE,discard0=FALSE,
                  ssRes=1,traplocs=traplocs[,-1],buffer=buffer,
                  encmod="P")
  simDat$S<-cbind(simDat$S,pop[,-1])
  
  ###2 - generate SC and SPIM data####
  #str(simDat)
  
  #2a - sc data is the counts of inds at a trap on a day: J (73) x K (92)
  #keep only the detected individuals from the simulated encounters
  # sc detections are summed across individuals, per site,
  
  dim(simDat$Y)
  max(simDat$Y)
  simSC<-simDat$Y
  simSC<-apply(simDat$Y,c(2,3),sum) 
  dim(simSC)
  max(simSC)
  
  hist_simSC$n[hist_simSC$Count%in%names(table(simSC))]<-
    hist_simSC$n[hist_simSC$Count%in%names(table(simSC))]+table(simSC)
  
  #2b - SPIM data
  #first, get the detections from the simulated encounters
  simSPIM_long<-as.data.frame(which(simDat$Y>0,arr.ind = TRUE))
  simSPIM_long$n<-NA
  for(i in 1:nrow(simSPIM_long)){
    simSPIM_long$n[i]<-simDat$Y[simSPIM_long[,1][i],simSPIM_long[,2][i],simSPIM_long[,3][i]]
  }
  simSPIM_long <- as.data.frame(lapply(simSPIM_long, rep, simSPIM_long$n))
  simSPIM_long<-simSPIM_long[,-ncol(simSPIM_long)]
  colnames(simSPIM_long)<-c("Ind","Trap","Occ")
  
  #add the partial IDs; the actual Ind column wont make its way into the final version
  #unique(simSPIM_long$Ind)
  simSPIM_long$Sex<-pop$Sex[match(simSPIM_long$Ind, pop$Ind)]
  simSPIM_long$Antlers<-pop$Antlers[match(simSPIM_long$Ind, pop$Ind)]
  simSPIM_long$Collar<-pop$Collar[match(simSPIM_long$Ind, pop$Ind)]
  
  #put various pieces into format for SPIM model
  #encounters
  y.obs<-array(0,dim=c(nrow(simSPIM_long),nrow(traplocs),K_days))
  for(ind in 1:nrow(simSPIM_long)){
    y.obs[ind,simSPIM_long$Trap[ind],simSPIM_long$Occ[ind]]<-1
  }
  #partial IDs, complete and partially missing (which are recorded as 0))
  #and for either using all 3 IDs or just sex and antlers.
  G.obs_complete_3ID<-cbind(simSPIM_long$Sex,simSPIM_long$Antler,simSPIM_long$Collar)
  G.obs_complete_2ID<-cbind(simSPIM_long$Sex,simSPIM_long$Antler)
  G.obs_missing_3ID<-G.obs_complete_3ID
  G.obs_missing_3ID[,1][sample(1:nrow(G.obs_complete_3ID),ceiling(nrow(G.obs_complete_3ID)/3))]<-0 #sex
  G.obs_missing_3ID[,2][sample(1:nrow(G.obs_complete_3ID),ceiling(nrow(G.obs_complete_3ID)/3))]<-0 #antlers
  G.obs_missing_2ID<-G.obs_missing_3ID[,-3]
  
  ncat_3ID<-ncol(G.obs_complete_3ID)
  ncat_2ID<-ncol(G.obs_complete_2ID)
  
  #Store unique genotypes that appear in the datasets
  IDcovs_complete_3ID<-vector("list",ncat_3ID)
  IDcovs_missing_3ID<-vector("list",ncat_3ID)
  IDcovs_complete_2ID<-vector("list",ncat_2ID)
  IDcovs_missing_2ID<-vector("list",ncat_2ID)
  
  for(c in 1:length(IDcovs_complete_3ID)){
    if(max(G.obs_complete_3ID[,c])==999){ #if working on antler cov and there are 999 values, ie with no antlers
      #because otherwise it will create categories all the way up to 999 which you dont want
      IDcovs_complete_3ID[[c]]<-c(1:(sort(unique(G.obs_complete_3ID[,c]), TRUE)[2]),999) 
    }else{
      IDcovs_complete_3ID[[c]]<-1:max(G.obs_complete_3ID[,c],na.rm = TRUE)  
    }
    
  }
  for(c in 1:length(IDcovs_missing_3ID)){
    if(max(G.obs_missing_3ID[,c])==999){ #if working on antler cov and there are 999 values, ie with no antlers
      #because otherwise it will create categories all the way up to 999 which you dont want
      IDcovs_missing_3ID[[c]]<-c(1:(sort(unique(G.obs_missing_3ID[,c]), TRUE)[2]),999) 
    }else{
      IDcovs_missing_3ID[[c]]<-1:max(G.obs_missing_3ID[,c],na.rm = TRUE)  
    }
    
  }
  for(c in 1:length(IDcovs_complete_2ID)){
    if(max(G.obs_complete_2ID[,c])==999){ #if working on antler cov and there are 999 values, ie with no antlers
      #because otherwise it will create categories all the way up to 999 which you dont want
      IDcovs_complete_2ID[[c]]<-c(1:(sort(unique(G.obs_complete_2ID[,c]), TRUE)[2]),999) 
    }else{
      IDcovs_complete_2ID[[c]]<-1:max(G.obs_complete_2ID[,c],na.rm = TRUE)  
    }
    
  }
  for(c in 1:length(IDcovs_missing_2ID)){
    if(max(G.obs_missing_2ID[,c])==999){ #if working on antler cov and there are 999 values, ie with no antlers
      #because otherwise it will create categories all the way up to 999 which you dont want
      IDcovs_missing_2ID[[c]]<-c(1:(sort(unique(G.obs_missing_2ID[,c]), TRUE)[2]),999) 
    }else{
      IDcovs_missing_2ID[[c]]<-1:max(G.obs_missing_2ID[,c],na.rm = TRUE)  
    }
  }
  
  IDlist_complete_3ID<-list(ncat=ncat_3ID,IDcovs=IDcovs_complete_3ID)
  IDlist_missing_3ID<-list(ncat=ncat_3ID,IDcovs=IDcovs_missing_3ID)
  IDlist_complete_2ID<-list(ncat=ncat_2ID,IDcovs=IDcovs_complete_2ID)
  IDlist_missing_2ID<-list(ncat=ncat_2ID,IDcovs=IDcovs_missing_2ID)
  
  #gamma, the initial values for pop level proportions for each catvar
  #based on their proportions in the data 
  gamma_props_priors_complete_3ID<-vector("list",ncat_3ID)
  gamma_props_priors_missing_3ID<-vector("list",ncat_3ID)
  gamma_props_priors_complete_2ID<-vector("list",ncat_2ID)
  gamma_props_priors_missing_2ID<-vector("list",ncat_2ID)
  
  for(c in 1:2){
    gamma_props_priors_complete_3ID[[c]]<-rep(1/length(IDcovs_complete_3ID[[c]]),length(IDcovs_complete_3ID[[c]]))
    gamma_props_priors_missing_3ID[[c]]<-rep(1/length(IDcovs_missing_3ID[[c]]),length(IDcovs_missing_3ID[[c]]))
    gamma_props_priors_complete_2ID[[c]]<-rep(1/length(IDcovs_complete_2ID[[c]]),length(IDcovs_complete_2ID[[c]]))
    gamma_props_priors_missing_2ID[[c]]<-rep(1/length(IDcovs_missing_2ID[[c]]),length(IDcovs_missing_2ID[[c]]))
  }
  gamma_props_priors_complete_3ID[[3]]<-c((100-length(which(G.obs_complete_3ID[,3]==2)))/100,length(which(G.obs_complete_3ID[,3]==2))/100)
  gamma_props_priors_missing_3ID[[3]]<-c((100-length(which(G.obs_missing_3ID[,3]==2)))/100,length(which(G.obs_missing_3ID[,3]==2))/100)
  
  #initial value for sigma
  sigma<-2.1
  
  #obstype  (poisson, bernoulli)
  obstype<-"poisson"
  
  # wrap data int a list
  dataSPIM_complete_3ID<-list(y.obs=y.obs,G.obs=G.obs_complete_3ID,IDlist=IDlist_complete_3ID,
                              X=traplocs[,-1],K=K_days,buff=buffer, obstype=obstype)
  dataSPIM_missing_3ID<-list(y.obs=y.obs,G.obs=G.obs_missing_3ID,IDlist=IDlist_missing_3ID,
                             X=traplocs,K=K_days, buff=buffer, obstype=obstype)
  dataSPIM_complete_2ID<-list(y.obs=y.obs,G.obs=G.obs_complete_2ID,IDlist=IDlist_complete_2ID,
                              X=traplocs[,-1],K=K_days,buff=buffer, obstype=obstype)
  dataSPIM_missing_2ID<-list(y.obs=y.obs,G.obs=G.obs_missing_2ID,IDlist=IDlist_missing_2ID,
                             X=traplocs,K=K_days, buff=buffer, obstype=obstype)
  
  print(it)
  print(paste0(nrow(simSPIM_long), " detections of ",
               length(unique(simSPIM_long$Ind)) ," inds"))
  
  datasets_list[[it]]<-list(traplocs=traplocs, buffer=buffer,
                            simDat=simDat,
                            simSC=simSC,
                            dataSPIM_complete_3ID=dataSPIM_complete_3ID,
                            dataSPIM_missing_3ID=dataSPIM_missing_3ID,
                            dataSPIM_complete_2ID=dataSPIM_complete_2ID,
                            dataSPIM_missing_2ID=dataSPIM_missing_2ID,
                            gamma_props_priors_complete_3ID=gamma_props_priors_complete_3ID,
                            gamma_props_priors_missing_3ID=gamma_props_priors_missing_3ID,
                            gamma_props_priors_complete_2ID=gamma_props_priors_complete_2ID,
                            gamma_props_priors_missing_2ID=gamma_props_priors_missing_2ID)
  
  
  hist_simSPIM <-c(hist_simSPIM ,nrow(simSPIM_long)) 
  hasCollars<-c(hasCollars,ifelse(length(dataSPIM_missing_3ID$IDlist[[2]][[3]])==2,"yes","no"))
}


#how do these simulations with this detection probability 
#compare with observed 2019 SPIM data
hist(hist_simSPIM,xlim=c(0,200))
abline(v=68,lty=2) #observed

#compared to breakdown of SC counts
hist_simSC$prop<-hist_simSC$n/sum(hist_simSC$n)
hist_simSC$Data<-"Simulated"
hist_obsSC<-data.frame(Count=c(0,1,2,4,8),n=c(6392,32,13,2,1))
hist_obsSC$prop<-hist_obsSC$n/sum(hist_obsSC)
hist_obsSC$Data<-"Observed"
hist_SC<-rbind(hist_simSC,hist_obsSC)
ggplot(hist_SC,aes(x = Count,y=prop,fill=Data))+
  geom_bar(stat="identity",position='dodge')

#keep just the 50 simulated datasets that most look like our data
#has the number of SPIM detections
#has collars
lower<-round(mean(hist_simSPIM)-(sd(hist_simSPIM)),0)
upper<-round(mean(hist_simSPIM)+(sd(hist_simSPIM)),0)
toKeep<-intersect(which(lower<hist_simSPIM& hist_simSPIM>upper),which(hasCollars=="yes"))
datasets_list<-datasets_list[toKeep]


save(datasets_list, file = "simDatasets_Poiss_forSCandSPIM.rda")