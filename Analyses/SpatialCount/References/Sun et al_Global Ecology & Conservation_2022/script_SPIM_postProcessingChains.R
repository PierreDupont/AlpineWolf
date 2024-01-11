#post processing SPIM model chains
library(coda)
library(tidyr)
library(ggplot2)
library(abind)
library(bayesplot)

#data and prep
study<-"Algar"
species1<-"caribou"
timeFrame<-"3mos"
dataType<-"all_wCollars"
year<- 2019

#some functions
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
getmode <- function(v) {
  uniqv <- unique(round(v,2))
  uniqv[which.max(tabulate(match(round(v,2), uniqv)))]
}


#bring each chain 
load(paste0("outputSPIM_",species1,"_",study,"_",year,"_",timeFrame,"_",dataType,"_chain1.RData"))
out_chain1<-out
#dataType<-"all_M800" #M400_wCollars" #all GrpSize1 all_wCollars
load(paste0("outputSPIM_",species1,"_",study,"_",year,"_",timeFrame,"_",dataType,"_chain2.RData"))
out_chain2<-out
dataType<-"all_wCollars" #M400_wCollars" #all GrpSize1 all_wCollars
load(paste0("outputSPIM_",species1,"_",study,"_",year,"_",timeFrame,"_",dataType,"_chain3.RData"))
out_chain3<-out
dataType<-"all_wCollars" #M400_wCollars" #all GrpSize1 all_wCollars

#how much to burn
nburn<-2000

#combine relevant parts of the output
out_chainsAll<-list(mcmc(out_chain1$out[-(1:nburn),]),
                    mcmc(out_chain2$out[-(1:nburn),]),
                    mcmc(out_chain3$out[-(1:nburn),]))
out_IDsAll<-abind(out_chain1$IDout[-(1:nburn),],out_chain2$IDout[-(1:nburn),],out_chain3$IDout[-(1:nburn),],
                  along=3)
out_chainsInOne<-rbind(out_chain1$out,out_chain2$out,out_chain3$out)

# rejection rates, convergence diagnostics, and effective sizes
1-rejectionRate(as.mcmc.list(out_chainsAll))
gel<-gelman.diag(out_chainsAll, autoburnin=FALSE,multivariate=FALSE)
neff<-effectiveSize(out_chainsAll)

gel
neff
assign(paste0("gel_",commonName,"_",year,"_",CapStr(dataType)),gel)
assign(paste0("neff_",commonName,"_",year,"_",CapStr(dataType)),neff)


#summaries
HPDI_95<-HPDinterval(as.mcmc(out_chainsInOne),prob=0.95)
HPDI_50<-HPDinterval(as.mcmc(out_chainsInOne),prob=0.50)
mode_estimates<-apply(out_chainsInOne,2,getmode)
summaryChains<-summary(as.mcmc.list(out_chainsAll))
HPDI_95
HPDI_50
mode_estimates
summaryChains
assign(paste0("HPDI95_",commonName,"_",year,"_",CapStr(dataType)),HPDI_95)
assign(paste0("HPDI50_",commonName,"_",year,"_",CapStr(dataType)),HPDI_50)
assign(paste0("mode_",commonName,"_",year,"_",CapStr(dataType)),mode_estimates)
assign(paste0("summaryChains_",species1,"_",study,"_",year,"_",dataType),summaryChains)

#plots
png(paste0("mcmc_",species1,"_Algar_",year,"_",timeFrame,"_",dataType,".png"),width = 720)
mcmc_hist(as.mcmc.list(out_chainsAll),pars=c("N","lam0","sigma","psi","n"))
dev.off()

niters=niter-nburn

PID<-matrix(0,nrow=dim(out_IDsAll)[2],ncol=1+dim(out_IDsAll)[2])
PID[,1]<-seq(1,dim(out_IDsAll)[2],1)
PID<-as.data.frame(PID)
colnames(PID)<-c("ID",seq(1,dim(out_IDsAll)[2],1))
rownames(PID)<-colnames(PID)[-1]

for(z in 1:ncol(out$IDout)){ #for each individual
  print(z)
  check=z #Which sample to check
  final_rowSums<-matrix(NA,ncol=dim(out_IDsAll)[2],nrow=3)
  
  for(c in 1:3){ #for each chain
    storematch=matrix(0,nrow=dim(out_IDsAll)[2],ncol=niters)
    
    for(i in 1:niters){ # for each iteration
      storematch[,i]=out_IDsAll[i,check,c]==out_IDsAll[i,,c]
    }# for an iteration in a chain
    tmp_rowSums<-rowSums(storematch) 
    final_rowSums[c,]<-tmp_rowSums
  }# for  chain
  
  
  PID[z,2:ncol(PID)]<-colSums(final_rowSums)/(niters*3)
}

PID<-PID%>%
  pivot_longer(-ID,names_to="ComparedTo_Ind",values_to="Probability")
PID$ComparedTo_Ind<-as.numeric(PID$ComparedTo_Ind)

plot_posteriorIDprob<-ggplot(PID, aes(x=ID, y=ComparedTo_Ind)) +
  geom_tile(aes(fill = Probability)) +
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title = paste0(study," ",year,": Posterior Pairwise Probability of Identity"),
       subtitle=paste0("Using data from ",dataType, " group sizes"),
       x="ID",y="ID")
plot_posteriorIDprob
ggsave(paste0("plotpostID_",species1,"_",study,"_",year,"_",dataType,".png"), plot = plot_posteriorIDprob)
assign(paste0("plotpostID_",species1,"_",study,"_",year,"_",dataType),plot_posteriorIDprob)

rm("out","out_chain1","out_chain2","out_chain3",
   "out_chainsAll", "out_chainsInOne","out_IDsAll",
   "gel","neff","HPDI_95","HPDI_50",
   "storematch","plot_posteriorIDprob")

save.image(file=paste0("outputSPIM_",LowStr(commonName),"_",study,"_",year,"_",dataType,"_Poiss_postProcessed.RData"))

