# Organizing Simulation results

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
library(dplyr)

#
#### SPIM simulation results#####
neff_all<-data.frame()
plot_data<-data.frame()

RData_SPIM<-list.files( pattern =  "_postProcessed.RData")

length(RData_SPIM)
tmp_list_SPIM<-list(names=RData_SPIM,
               plotN=plot_data,plotn=plot_data,
               plotLam0=plot_data,plotSigma=plot_data,
               neff=neff_all)

#to get mean median mode and hpdi instead of using bayesplot
for(r in 1:length(RData_SPIM)){ #
  print(RData_SPIM[r])
  load(RData_SPIM[r])
  
  HPDI95<-HPDI_95
  HPDI50<-HPDI_50
  summaryChains<-summaryChains
  
  
  posteriorValues_N<-c(HPDI95[3,1],HPDI50[3,1],HPDI50[3,2],HPDI95[3,2],
                       summaryChains$statistics[,1][3],summaryChains$quantiles[,3][3],
                       mode_estimates[3],
                       summaryChains$statistics[,2][which(rownames(summaryChains$statistics)=="N")])
  posteriorValues_n<-c(HPDI95[4,1],HPDI50[4,1],HPDI50[4,2],HPDI95[4,2],
                       summaryChains$statistics[,1][4],summaryChains$quantiles[,3][4],
                       mode_estimates[4],
                       summaryChains$statistics[,2][which(rownames(summaryChains$statistics)=="n")])
  posteriorValues_Lam0<-c(HPDI95[1,1],HPDI50[1,1],HPDI50[1,2],HPDI95[1,2],
                          summaryChains$statistics[,1][1],summaryChains$quantiles[,3][1],
                          mode_estimates[1],
                          summaryChains$statistics[,2][which(rownames(summaryChains$statistics)=="lam0")])
  posteriorValues_Sigma<-c(HPDI95[2,1],HPDI50[2,1],HPDI50[2,2],HPDI95[2,2],
                           summaryChains$statistics[,1][2],summaryChains$quantiles[,3][2],
                           mode_estimates[2],
                           summaryChains$statistics[,2][which(rownames(summaryChains$statistics)=="sigma")])
  
  neff_tmp<-neff 
  tmp_neff<-t(as.data.frame(neff_tmp))
  
  tmp_list_SPIM$plotN<-rbind(tmp_list_SPIM$plotN,posteriorValues_N)
  tmp_list_SPIM$plotn<-rbind(tmp_list_SPIM$plotn,posteriorValues_n)
  tmp_list_SPIM$plotLam0<-rbind(tmp_list_SPIM$plotLam0,posteriorValues_Lam0) 
  tmp_list_SPIM$plotSigma<-rbind(tmp_list_SPIM$plotSigma,posteriorValues_Sigma)
  
  tmp_list_SPIM$neff<-rbind(tmp_list_SPIM$neff,tmp_neff)
  tmp_list_SPIM$names[r]<-RData_SPIM[r]
  
}

colnames(tmp_list_SPIM$plotN)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list_SPIM$plotn)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list_SPIM$plotLam0)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list_SPIM$plotSigma)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")

tmp_list_SPIM$plotN$CoV<-tmp_list_SPIM$plotN$SD/tmp_list_SPIM$plotN$Mean
tmp_list_SPIM$plotn$CoV<-tmp_list_SPIM$plotn$SD/tmp_list_SPIM$plotn$Mean
tmp_list_SPIM$plotLam0$CoV<-tmp_list_SPIM$plotLam0$SD/tmp_list_SPIM$plotLam0$Mean
tmp_list_SPIM$plotSigma$CoV<-tmp_list_SPIM$plotSigma$SD/tmp_list_SPIM$plotSigma$Mean

tmp_list_SPIM$plotN$Iter<-as.numeric(sapply(strsplit(RData_SPIM,"_"),"[[",4))
tmp_list_SPIM$plotLam0$Iter<-as.numeric(sapply(strsplit(RData_SPIM,"_"),"[[",4))
tmp_list_SPIM$plotSigma$Iter<-as.numeric(sapply(strsplit(RData_SPIM,"_"),"[[",4))


#### SC SIM RESULTS######

getmode <- function(v) {
  uniqv <- unique(round(v,2))
  uniqv[which.max(tabulate(match(round(v,2), uniqv)))]
}

neff_all<-data.frame()
mcmc_plot<-data.frame()
plot_data<-data.frame()

RData_SC<-list.files( pattern =  "outputNimble_simSC_")

length(RData_SC)
tmp_list_SC<-list(names=RData_SC,
               plotN=plot_data,plotLam0=plot_data,plotSigma=plot_data,
               neff=neff_all)

#bring in sc results
for(r in 1:length(RData_SC)){
  print(RData_SC[r])
  load(RData_SC[r])
  if(exists("out.mcmc.bind_final")){
    tmp_out<-out.mcmc.bind_final  
  } else{
    tmp_out<-out.mcmc
  }
  
  tmp_out<-tmp_out[-seq(1,5000,1),]
  
  png(paste0("mcmc_simSC_it",r,".png"))
  print(mcmc_hist(tmp_out, pars = c("N","lam0","sigma","psi")))
  dev.off()
  
  
  tmp_out_inOne<-as.mcmc(rbind(tmp_out))
  hpdint_95<-HPDinterval(tmp_out_inOne,prob=0.95)
  hpdint_50<-HPDinterval(tmp_out_inOne,prob=0.50)
  
  whichD<-which(rownames(hpdint_50)=="D")
  whichD2<-which(colnames(tmp_out_inOne)=="D")
  
  whichN<-which(rownames(hpdint_50)=="N")
  whichN2<-which(colnames(tmp_out_inOne)=="N")
  
  whichlam0<-which(rownames(hpdint_50)=="lam0")
  whichlam02<-which(colnames(tmp_out_inOne)=="lam0")
  
  whichPsi<-which(rownames(hpdint_50)=="psi")
  whichPsi<-which(colnames(tmp_out_inOne)=="psi")
  
  whichSigma<-which(rownames(hpdint_50)=="sigma")
  whichSigma2<-which(colnames(tmp_out_inOne)=="sigma")
  
  posteriorValues_N<-c(hpdint_95[whichN,1],hpdint_50[whichN,1],hpdint_50[whichN,2],hpdint_95[whichN,2],
                       mean(tmp_out_inOne[,whichN2]),median(tmp_out_inOne[,whichN2]),getmode(tmp_out_inOne[,whichN2]),
                       sd(tmp_out_inOne[,whichN2]))
  posteriorValues_lam0<-c(hpdint_95[whichlam0,1],hpdint_50[whichlam0,1],hpdint_50[whichlam0,2],hpdint_95[whichlam0,2],
                          mean(tmp_out_inOne[,whichlam02]),median(tmp_out_inOne[,whichlam02]),getmode(tmp_out_inOne[,whichlam02]),
                          sd(tmp_out_inOne[,whichlam02]))
  posteriorValues_Sigma<-c(hpdint_95[whichSigma,1],hpdint_50[whichSigma,1],hpdint_50[whichSigma,2],hpdint_95[whichSigma,2],
                           mean(tmp_out_inOne[,whichSigma2]),median(tmp_out_inOne[,whichSigma2]),getmode(tmp_out_inOne[,whichSigma2]),
                           sd(tmp_out_inOne[,whichSigma2]))
  
  tmp_neff<-t(as.data.frame(neff[c(whichD,whichN,whichlam0,whichPsi,whichSigma)]))
  
  tmp_list_SC$plotN<-rbind(tmp_list_SC$plotN,posteriorValues_N)
  tmp_list_SC$plotLam0<-rbind(tmp_list_SC$plotLam0,posteriorValues_lam0)
  tmp_list_SC$plotSigma<-rbind(tmp_list_SC$plotSigma,posteriorValues_Sigma)
  
  tmp_list_SC$neff<-rbind(tmp_list_SC$neff,tmp_neff)
}

colnames(tmp_list_SC$plotN)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list_SC$plotLam0)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list_SC$plotSigma)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")

tmp_list_SC$plotN$CoV<-tmp_list_SC$plotN$SD/tmp_list_SC$plotN$Mean
tmp_list_SC$plotLam0$CoV<-tmp_list_SC$plotLam0$SD/tmp_list_SC$plotLam0$Mean
tmp_list_SC$plotSigma$CoV<-tmp_list_SC$plotSigma$SD/tmp_list_SC$plotSigma$Mean

tmp_list_SC$plotN$Iter<-as.numeric(sapply(strsplit(sapply(strsplit(RData_SC,"_"),"[[",3),".R"),"[[",1))
tmp_list_SC$plotLam0$Iter<-as.numeric(sapply(strsplit(sapply(strsplit(RData_SC,"_"),"[[",3),".R"),"[[",1))
tmp_list_SC$plotSigma$Iter<-as.numeric(sapply(strsplit(sapply(strsplit(RData_SC,"_"),"[[",3),".R"),"[[",1))

rm(list=setdiff(ls(), c("data_perIt","RData_SPIM","tmp_list_SPIM","r","inDir2","getmode",
                        "RData_SC","tmp_list_SC")))


### prep for plots ####

N_long_SPIM<- reshape2::melt(tmp_list_SPIM$plotN, id.vars = c("Iter"),measure.vars = c(1:7,9))
N_long_SPIM$Model<-"SPIM"

N_long_SC<- reshape2::melt(tmp_list_SC$plotN, id.vars = c("Iter"), measure.vars = c(1:7,9))
N_long_SC$Model<-"SC"

if(sum(colnames(N_long_SC)==colnames(N_long_SPIM))==ncol(N_long_SPIM)){
  N_long<-rbind(N_long_SC,N_long_SPIM)  
}


Sigma_long_SPIM<- reshape2::melt(tmp_list_SPIM$plotSigma, id.vars = c("Iter"),measure.vars = c(1:7,9))
Sigma_long_SPIM$Model<-"SPIM"

Sigma_long_SC<- reshape2::melt(tmp_list_SC$plotSigma, id.vars = c("Iter"), measure.vars = c(1:7,9))
Sigma_long_SC$Model<-"SC"

if(sum(colnames(Sigma_long_SC)==colnames(Sigma_long_SPIM))==ncol(Sigma_long_SPIM)){
  Sigma_long<-rbind(Sigma_long_SC,Sigma_long_SPIM)  
}

# spim plots
plotN_SPIM<-ggplot(data=N_long_SPIM[N_long_SPIM$variable%in%c("Median","Mean","Mode"),])+
  geom_point(aes(x=variable,y=value))+
  geom_hline(yintercept=141,lty=2)+
  labs(title="SPIM Abundance with different point estimates")
plotN_SPIM

plotN_SPIM<-ggplot(N_long_SPIM[N_long_SPIM$variable%in%c("Median","Mean","Mode"),], aes(x=variable, y=value)) + 
  geom_boxplot()+
  geom_hline(yintercept=141,lty=2)+
  labs(title="SPIM Abundance with different point estimates")
plotN_SPIM


#sc plots
plotN_SC<-ggplot(data=N_long_SC[N_long_SC$variable%in%c("Median","Mean","Mode"),])+
  geom_point(aes(x=variable,y=value))+
  geom_hline(yintercept=141,lty=2)+
  labs(title="SC Abundance with different point estimates")
plotN_SC

plotN_SC<-ggplot(N_long_SC[N_long_SC$variable%in%c("Median","Mean","Mode"),], aes(x=variable, y=value)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  geom_hline(yintercept=141,lty=2)+
  labs(title="SC Abundance with different point estimates")
plotN_SC


#all together
plotN<-ggplot(N_long[N_long$variable%in%c("Median","Mean","Mode"),], 
              aes(x=variable, y=value,color=Model)) + 
  geom_boxplot()+
  geom_hline(yintercept=141,lty=2)+
  labs(title="Spatial Count v. Spatial Partial Identity Models",
       x="Point Estimate",y="Estimated Abundance")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plotN

plotSigma<-ggplot(Sigma_long[Sigma_long$variable%in%c("Median","Mean","Mode"),], 
       aes(x=variable, y=value,color=Model)) + 
  geom_boxplot()+
  geom_hline(yintercept=2.8,lty=2)+
  labs(title="Spatial Count v. Spatial Partial Identity Models",
       x="Point Estimate",y="Estimated Sigma")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


####Bias ####
N_long_ptEst<-N_long[N_long$variable %in% c("Mean","Median","Mode","CoV"),]
N_long_ptEst$variable<-factor(N_long_ptEst$variable,levels = c("Mode","Median","Mean","CoV"))
N_long_ptEst$Bias<-(N_long_ptEst$value-137)/137
N_long_ptEst$Param<-"N"


n<-Sigma_long[Sigma_long$variable %in% c("Mean","Median","Mode","CoV"),]
Sigma_long_ptEst$variable<-factor(Sigma_long_ptEst$variable,levels = c("Mode","Median","Mean","CoV"))
Sigma_long_ptEst$Bias<-(Sigma_long_ptEst$value-2.1)/2.1 #2.1
Sigma_long_ptEst$Param<-"Sigma"
Bias_long<-rbind(N_long_ptEst,Sigma_long_ptEst)
Bias_long$Bias<-Bias_long$Bias*100


plot_Bias_modeOnly<-ggplot(Bias_long[Bias_long$variable=="Mode",], 
                  aes(x=Param, y=Bias,fill=Model)) + 
  geom_boxplot()+
  geom_hline(yintercept=0,lty=2)+
  labs(x="Parameter",y="Bias (%)")+
  scale_fill_manual(values=c("blue","red"))+#c( "#0072B2","#999999"))+
  theme(legend.position = "bottom",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA)
        )
plot_Bias_modeOnly

###Coverage ####
#Abundance
for(i in 1:nrow(tmp_list_SPIM$plotN)){
  tmp_list_SPIM$plotN$Coverage[i]<-ifelse(tmp_list_SPIM$plotN$HPDI_5[i]<137&tmp_list_SPIM$plotN$HPDI_95[i]>137,1,0)
}
sum(tmp_list_SPIM$plotN$Coverage)/nrow(tmp_list_SPIM$plotN)

for(i in 1:nrow(tmp_list_SC$plotN)){
  tmp_list_SC$plotN$Coverage[i]<-ifelse(tmp_list_SC$plotN$HPDI_5[i]<137&tmp_list_SC$plotN$HPDI_95[i]>137,1,0)
}
sum(tmp_list_SC$plotN$Coverage)/nrow(tmp_list_SC$plotN)

#sigma
for(i in 1:nrow(tmp_list_SPIM$plotSigma)){
  tmp_list_SPIM$plotSigma$Coverage[i]<-ifelse(tmp_list_SPIM$plotSigma$HPDI_5[i]<2.1&tmp_list_SPIM$plotSigma$HPDI_95[i]>2.1,1,0)
}
sum(tmp_list_SPIM$plotSigma$Coverage)/nrow(tmp_list_SPIM$plotSigma)

for(i in 1:nrow(tmp_list_SC$plotSigma)){
  tmp_list_SC$plotSigma$Coverage[i]<-ifelse(tmp_list_SC$plotSigma$HPDI_5[i]<2.1&tmp_list_SC$plotSigma$HPDI_95[i]>2.1,1,0)
}
sum(tmp_list_SC$plotSigma$Coverage)/nrow(tmp_list_SC$plotSigma)



####Coefficient of Variation ####

plot_CoV_N<-ggplot(N_long_ptEst[N_long_ptEst$variable=="CoV",], 
       aes(x = Model, y = value*100,label=Iter)) +
  # geom_jitter(position = position_jitter(width = 0.1, height = 0))+
  geom_line(aes(group = Iter),lwd=1.5,alpha=0.5)+#,lty=Year
  geom_point(size=4,pch=20, color="black")+
  # geom_text()+
#  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 20))+#breaks = seq(0, 180, by = 20))+
 # facet_wrap(vars(Study))+
  labs(x="",y="Coefficient of Variation (%)")+
  theme(legend.position = "right",
        legend.justification = "center",
        legend.key=element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x=element_blank(),
        plot.title = element_text(size=20, hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))




plot_CoV_Sigma<-ggplot(Sigma_long_ptEst[Sigma_long_ptEst$variable=="CoV",], 
       aes(x = Model, y = value*100,label=Iter)) +
  # geom_jitter(position = position_jitter(width = 0.1, height = 0))+
  geom_line(aes(group = Iter),lwd=1.5,alpha=0.5)+#,lty=Year
  geom_point(size=4,pch=20, color="black")+
  # geom_text()+
  #  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 20))+#breaks = seq(0, 180, by = 20))+
  # facet_wrap(vars(Study))+
  labs(x="",y="Coefficient of Variation (%)")+
  theme(legend.position = "right",
        legend.justification = "center",
        legend.key=element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x=element_blank(),
        plot.title = element_text(size=20, hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))

#save simulated COVs to add to the empirical results
tmp<-rbind(N_long_ptEst,Sigma_long_ptEst)
tmp%>%group_by(Model,Param,variable)%>%summarize(mean=mean(value),
                                             min=min(value),
                                             max=max(value))

sim_long<-rbind(N_long_ptEst[N_long_ptEst$variable=="CoV",],
      Sigma_long_ptEst[Sigma_long_ptEst$variable=="CoV",])
sim_long$Bias[which(sim_long$variable=="CoV")]<-NA
sim_long%>%group_by(Model,Param)%>%summarize(mean=mean(value),
                                       min=min(value),
                                       max=max(value))
save(sim_long, file = "simulated_CoV.rda")




plot_fig4<-plot_Bias_modeOnly
