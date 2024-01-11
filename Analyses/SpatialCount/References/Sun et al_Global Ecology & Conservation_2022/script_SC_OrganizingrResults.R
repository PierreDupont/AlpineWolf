#
#organize SC results
#

### --- Load packages, prep setting ####
library(data.table)
library(R2jags)
library(coda)
library(bayesplot)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(scales)
library(lubridate)
library(stringr)
library(ggrepel)

getmode <- function(v) {
  uniqv <- unique(round(v,2))
  uniqv[which.max(tabulate(match(round(v,2), uniqv)))]
}

# SC dataprep ####
neff_all<-data.frame()
rhat_all<-data.frame()
rhat_allUpper<-data.frame()
indsDetectedMaxPossible_final<-c()
indsDetectedPerDay_max_final<-c()
daySitesWithDets_final<-data.frame()
plot_list<-list()
mcmc_plot<-data.frame()
plot_data<-data.frame()

RData<-c("outputNimble_caribou_Algar_2016_3mos_sum.RData",
         "outputNimble_caribou_Algar_2017_3mos_sum.RData",
         "outputNimble_caribou_Algar_2018_3mos_sum.RData",
         "outputNimble_caribou_Algar_2019_3mos_sum.RData")

length(RData)
tmp_list<-list(names=RData,
               plotN=plot_data,plotLam0=plot_data,plotSigma=plot_data,
               indsDetectedMaxPossible=indsDetectedMaxPossible_final,
               indsDetectedPerDay=indsDetectedPerDay_max_final,
               daySitesWithDets=daySitesWithDets_final,
               neff=neff_all,rhat=rhat_all,rhat_Upper=rhat_allUpper)
#to get mean median mode and hpdi instead of using bayesplot
for(r in 1:length(RData)){ #
  print(RData[r])
  load(RData[r])
  if(exists("out.mcmc.bind_final")){
    tmp_out<-out.mcmc.bind_final  
  } else{
    tmp_out<-out.mcmc.bind
  }
  
  for(i in 1:length(tmp_out)){
    tmp_out[[i]]<-tmp_out[[i]][-seq(1,1000,1),]
  }
  
  png(paste0("mcmc_",substring(RData[r],14,nchar(RData[r])-6),".png"))
  print(mcmc_hist(tmp_out, pars = c("N","lam0","sigma","psi")))
  dev.off()

  #subset posteriors draws to the focal parameters
  params<-which(colnames(tmp_out[[1]])%in% params[1:5])
  tmp_out[[1]]<-tmp_out[[1]][,params]
  tmp_out[[2]]<-tmp_out[[2]][,params]
  tmp_out[[3]]<-tmp_out[[3]][,params]
  
  #get summaries
  summaryChains<-summary(mcmc.list(as.mcmc(tmp_out[[1]]),as.mcmc(tmp_out[[2]]),as.mcmc(tmp_out[[3]])))
  
  tmp_out_inOne<-as.mcmc(rbind(tmp_out[[1]],tmp_out[[2]],tmp_out[[3]]))
  hpdint_95<-HPDinterval(tmp_out_inOne,prob=0.95)
  hpdint_50<-HPDinterval(tmp_out_inOne,prob=0.50)
  
  posteriorValues_N<-c(hpdint_95[2,1],hpdint_50[2,1],hpdint_50[2,2],hpdint_95[2,2],
                       mean(tmp_out_inOne[,2]),median(tmp_out_inOne[,2]),getmode(tmp_out_inOne[,2]),
                       summaryChains$statistics[,2][which(rownames(summaryChains$statistics)=="N")])
  
  posteriorValues_lam0<-c(hpdint_95[3,1],hpdint_50[3,1],hpdint_50[3,2],hpdint_95[3,2],
                          mean(tmp_out_inOne[,3]),median(tmp_out_inOne[,3]),getmode(tmp_out_inOne[,3]),
                          summaryChains$statistics[,2][which(rownames(summaryChains$statistics)=="lam0")])
  
  posteriorValues_Sigma<-c(hpdint_95[which(rownames(hpdint_95)=="sigma"),1],hpdint_50[which(rownames(hpdint_95)=="sigma"),1],
                           hpdint_50[which(rownames(hpdint_95)=="sigma"),2],hpdint_95[which(rownames(hpdint_95)=="sigma"),2],
                           mean(tmp_out_inOne[,which(rownames(hpdint_50)=="sigma")]),
                           median(tmp_out_inOne[,which(rownames(hpdint_50)=="sigma")]),
                           getmode(tmp_out_inOne[,which(rownames(hpdint_50)=="sigma")]),
                           summaryChains$statistics[,2][which(rownames(summaryChains$statistics)=="sigma")])
  
  tmp_indsDetectedMaxPossible<-sum(SC_data,na.rm = TRUE)
  tmp_indsDetectedPerDay<-max(SC_data,na.rm = TRUE)
  tmp_daySitesWithDets<- c(length(which(SC_data>0)),
                           length(which(colSums(SC_data)>0)),
                           length(which(rowSums(SC_data)>0)))
 
  if(exists("neff")){
    tmp_neff<-t(as.data.frame(neff[c(1,2,3,4,which(names(neff)=="sigma"))]))
    tmp_rhat<-gel$psrf[c(1,2,3,4,which(names(neff)=="sigma")),]
  } else{
    gel<-gelman.diag(tmp_out, autoburnin=FALSE,multivariate=FALSE)
    neff<-effectiveSize(tmp_out)
    
  }
  tmp_neff<-t(as.data.frame(neff[c(1,2,3,4,which(names(neff)=="sigma"))]))
  tmp_rhat<-gel$psrf[c(1,2,3,4,which(names(neff)=="sigma")),]
  
  tmp_list$plotN<-rbind(tmp_list$plotN,posteriorValues_N)
  tmp_list$plotLam0<-rbind(tmp_list$plotLam0,posteriorValues_lam0)
  tmp_list$plotSigma<-rbind(tmp_list$plotSigma,posteriorValues_Sigma)
  
  tmp_list$indsDetectedMaxPossible[r]<-tmp_indsDetectedMaxPossible
  tmp_list$indsDetectedPerDay[r]<-tmp_indsDetectedPerDay
  tmp_list$daySitesWithDets<-rbind(tmp_list$daySitesWithDets,tmp_daySitesWithDets)
  tmp_list$neff<-rbind(tmp_list$neff,tmp_neff)
  tmp_list$rhat<-rbind(tmp_list$rhat,t(tmp_rhat[,1]))
  tmp_list$rhat_Upper<-rbind(tmp_list$rhat_Upper,t(tmp_rhat[,2]))
  
  rm(list=ls(pattern="^out"))
  #
  
}
colnames(tmp_list$plotN)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD" )
colnames(tmp_list$plotLam0)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list$plotSigma)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list$daySitesWithDets)<-c("DaySites","Days","Sites")

tmp_list$plotN$CoV<-tmp_list$plotN$SD/tmp_list$plotN$Mean
tmp_list$plotLam0$CoV<-tmp_list$plotLam0$SD/tmp_list$plotLam0$Mean
tmp_list$plotSigma$CoV<-tmp_list$plotSigma$SD/tmp_list$plotSigma$Mean

save(tmp_list, file = "output_SCresults_informativePriors.RData")

