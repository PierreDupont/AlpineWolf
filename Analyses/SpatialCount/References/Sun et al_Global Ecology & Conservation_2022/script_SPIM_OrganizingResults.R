#
#organize SPIM results
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
# SPIM dataprep ####
neff_all<-data.frame()
rhat_all<-data.frame()
rhat_allUpper<-data.frame()
plot_data<-data.frame()


RData<-c("outputSPIM_caribou_Algar_2016_grpSize1_Poiss_postProcessed.RData",
         "outputSPIM_caribou_Algar_2017_GrpSize1_Poiss_postProcessed.RData",
         "outputSPIM_caribou_Algar_2018_GrpSize1_Poiss_postProcessed.RData",
         "outputSPIM_caribou_Algar_2019_grpSize1_Poiss_postProcessed.RData",
         "outputSPIM_caribou_Algar_2016_all_Poiss_postProcessed.RData",
         "outputSPIM_caribou_Algar_2017_all_Poiss_postProcessed.RData",
         "outputSPIM_caribou_Algar_2018_all_Poiss_postProcessed.RData",
         "outputSPIM_caribou_Algar_2019_all_wCollars_Poiss_postProcessed.RData")

length(RData)
tmp_list<-list(names=RData,
               plotN=plot_data,plotn=plot_data,
               plotLam0=plot_data,plotSigma=plot_data,
               neff=neff_all,rhat=rhat_all,rhat_Upper=rhat_allUpper)

#to get mean median mode and hpdi instead of using bayesplot
for(r in 1:length(RData)){ #
  print(RData[r])
  load(RData[r])
  
 # dataType<-"all"
  
  HPDI95<-get(paste0("HPDI95_",commonName,"_",year,"_",CapStr(dataType)))
  HPDI50<-get(paste0("HPDI50_",commonName,"_",year,"_",CapStr(dataType)))
  
  summaryChains<-get(paste0("summaryChains_",species1,"_",study,"_",year,"_",dataType))
 
   posteriorValues_N<-c(HPDI95[rownames(HPDI95)=="N",1],
                        HPDI50[rownames(HPDI50)=="N",1],
                        HPDI50[rownames(HPDI50)=="N",2],
                        HPDI95[rownames(HPDI95)=="N",2],
                       summaryChains$statistics[,1][rownames(summaryChains$statistics)=="N"],
                       summaryChains$quantiles[,3][rownames(summaryChains$quantiles)=="N"],
                       mode_estimates[attributes(mode_estimates)$name=="N"],
                       summaryChains$statistics[,2][rownames(summaryChains$statistics)=="N"])
   
  posteriorValues_n<-c(HPDI95[rownames(HPDI95)=="n",1],HPDI50[rownames(HPDI50)=="n",1],
                       HPDI50[rownames(HPDI50)=="n",2],HPDI95[rownames(HPDI95)=="n",2],
                       summaryChains$statistics[,1][rownames(summaryChains$statistics)=="n"],
                       summaryChains$quantiles[,3][rownames(summaryChains$statistics)=="n"],
                       mode_estimates[attributes(mode_estimates)$names=="n"],
                       summaryChains$statistics[,2][rownames(summaryChains$statistics)=="n"])
  
  posteriorValues_Lam0<-c(HPDI95[rownames(HPDI95)=="lam0",1],HPDI50[rownames(HPDI50)=="lam0",1],
                          HPDI50[rownames(HPDI50)=="lam0",2],HPDI95[rownames(HPDI95)=="lam0",2],
                       summaryChains$statistics[,1][rownames(summaryChains$statistics)=="lam0"],
                       summaryChains$quantiles[,3][rownames(summaryChains$statistics)=="lam0"],
                       mode_estimates[attributes(mode_estimates)$names=="lam0"],
                       summaryChains$statistics[,2][rownames(summaryChains$statistics)=="lam0"])
  
  posteriorValues_Sigma<-c(HPDI95[rownames(HPDI95)=="sigma",1],HPDI50[rownames(HPDI50)=="sigma",1],
                           HPDI50[rownames(HPDI50)=="sigma",2],HPDI95[rownames(HPDI95)=="sigma",2],
                           summaryChains$statistics[,1][rownames(summaryChains$statistics)=="sigma"],
                           summaryChains$quantiles[,3][rownames(summaryChains$statistics)=="sigma"],
                           mode_estimates[attributes(mode_estimates)$names=="sigma"],
                           summaryChains$statistics[,2][rownames(summaryChains$statistics)=="sigma"])
  
  if(dataType=="All_M800"){
    dataType<-"all"
  }
  neff_tmp<-get(paste0("neff_",commonName,"_",year,"_",CapStr(dataType)))
  rhat_tmp<-get(paste0("gel_",commonName,"_",year,"_",CapStr(dataType)))
  tmp_neff<-t(as.data.frame(neff_tmp))
  tmp_rhat<-rhat_tmp$psrf
  
  tmp_list$plotN<-rbind(tmp_list$plotN,posteriorValues_N)
  tmp_list$plotn<-rbind(tmp_list$plotn,posteriorValues_n)
  tmp_list$plotLam0<-rbind(tmp_list$ plotLam0,posteriorValues_Lam0) 
  tmp_list$plotSigma<-rbind(tmp_list$plotSigma,posteriorValues_Sigma)
  
  tmp_list$neff<-rbind(tmp_list$neff,tmp_neff)
  tmp_list$rhat<-rbind(tmp_list$rhat,t(tmp_rhat[,1]))
  tmp_list$rhat_Upper<-rbind(tmp_list$rhat_Upper,t(tmp_rhat[,2]))
  tmp_list$names[r]<-RData[r]
  
  rm(list=setdiff(ls(), c("RData","r","tmp_list","inDir2")))
}


colnames(tmp_list$plotN)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list$plotn)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list$plotLam0)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")
colnames(tmp_list$plotSigma)<-c("HPDI_5","HPDI_25","HPDI_50","HPDI_95","Mean","Median","Mode","SD")


tmp_list$plotN$CoV<-tmp_list$plotN$SD/tmp_list$plotN$Mean
tmp_list$plotn$CoV<-tmp_list$plotn$SD/tmp_list$plotn$Mean
tmp_list$plotLam0$CoV<-tmp_list$plotLam0$SD/tmp_list$plotLam0$Mean
tmp_list$plotSigma$CoV<-tmp_list$plotSigma$SD/tmp_list$plotSigma$Mean

save(tmp_list, file = "output_SPIMresults_Informative.RData")

### SPIM: create density plots with Cats results####
names_split<-str_split_fixed(tmp_list$names, "_",n=8)

tmp_list$plotN$Year<-substring(tmp_list$names,26,29)
tmp_list$plotN$`DataType`<-names_split[,5]
#tmp_list$plotN$`DataType`[which(names_split[,7]=="wCollars")]<-"all_withCollar"
tmp_list$plotN$`DataType`<-unlist(lapply(tmp_list$plotN$`DataType`,LowStr))
tmp_list$plotN$`DataType`<-factor(tmp_list$plotN$`DataType`,levels=c("grpSize1","all","all_withCollar"))
tmp_list$plotN


tmp_list$plotn$Year<-substring(tmp_list$names,26,29)
tmp_list$plotn$`DataType`<-names_split[,5]
#tmp_list$plotn$`DataType`[which(names_split[,7]=="wCollars")]<-"all_withCollar"
tmp_list$plotn$`DataType`<-unlist(lapply(tmp_list$plotn$`DataType`,LowStr))
tmp_list$plotn$`DataType`<-factor(tmp_list$plotn$`DataType`,levels=c("grpSize1","all","all_withCollar"))
tmp_list$plotn

tmp_list$plotSigma$Year<-substring(tmp_list$names,26,29)
tmp_list$plotSigma$`DataType`<-names_split[,5]
#tmp_list$plotSigma$`DataType`[which(names_split[,7]=="wCollars")]<-"all_withCollar"
tmp_list$plotSigma$`DataType`<-unlist(lapply(tmp_list$plotSigma$`DataType`,LowStr))
tmp_list$plotSigma$`DataType`<-factor(tmp_list$plotSigma$`DataType`,levels=c("grpSize1","all","all_withCollar"))
tmp_list$plotSigma

tmp_list$plotLam0$Year<-substring(tmp_list$names,26,29)
tmp_list$plotLam0$`DataType`<-names_split[,5]
#tmp_list$plotLam0$`DataType`[which(names_split[,7]=="wCollars")]<-"all_withCollar"
tmp_list$plotLam0$`DataType`<-unlist(lapply(tmp_list$plotLam0$`DataType`,LowStr))
tmp_list$plotLam0$`DataType`<-factor(tmp_list$plotLam0$`DataType`,levels=c("grpSize1","all","all_withCollar"))
tmp_list$plotLam0


#posterior estimates 
#abundance
plotN_withCat<-ggplot(data=tmp_list$plotN) + 
  geom_segment(  aes(x=DataType,y=HPDI_5,xend=DataType, yend=HPDI_95), size=1, color="lightblue")+
  geom_segment(  aes(x=DataType, y=HPDI_25,xend=DataType, yend=HPDI_50), size=2, color="darkblue")+
  geom_point( aes(x=DataType, y=Mode), size=3, shape=21, fill="white") +
  #geom_point( aes(x=DataType, y=Mean), size=3, shape=21, fill="red") +
  geom_text_repel(nudge_x=0.5,aes(x = DataType,  y = Mode,  label = round(Mode,digits = 0)))+
  coord_cartesian(ylim=c(0, 800))+ 
  facet_grid(cols=vars(Year),scales = "free_x")+ #
  theme(plot.title = element_text(hjust = 0.5))+#,strip.placement = "outside")+
  labs(title="Caribou Abundance in Algar using Spatial Partial ID Models",
       x ="Data type", y = "Abundance estimate")
plotN_withCat

# Number of observed inds
plotn_withCat<-ggplot(data=tmp_list$plotn) + 
  geom_segment(  aes(x=DataType,y=HPDI_5,xend=DataType, yend=HPDI_95), size=1, color="lightblue")+
  geom_segment(  aes(x=DataType, y=HPDI_25,xend=DataType, yend=HPDI_50), size=2, color="darkblue")+
  geom_point( aes(x=DataType, y=Mode), size=3, shape=21, fill="white") +
  geom_text_repel(nudge_x=0.5,aes(x = DataType,  y = Mode,  label = round(Mode,digits = 0)))+
  coord_cartesian(ylim=c(0, 50))+ 
  facet_grid(cols=vars(Year),scales = "free_x")+ #
  theme(plot.title = element_text(hjust = 0.5))+#,strip.placement = "outside")+
  labs(title="Estimated number of observed individuals using Spatial Partial ID Models",
       x ="Data type", y = "Observed Individuals")
plotn_withCat

SPIM_results<-tmp_list

save.image(file = "output_SPIMresults_Informative.RData")
