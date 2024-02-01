#combining/comparing SC and SPIM results for caribou
library(ggplot2)
library(ggrepel)
library(stringr)
library(rphylopic)
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

#####Load in Informed SC results ####
load("output_SCresults_informativePriors.RData")
SC_results_informed<-tmp_list
SC_results_informed$plotLam0$CoV<-SC_results_informed$plotLam0$SD/SC_results_informed$plotLam0$Mean
SC_results_informed$plotN$CoV<-SC_results_informed$plotN$SD/SC_results_informed$plotN$Mean
SC_results_informed$plotSigma$CoV<-SC_results_informed$plotSigma$SD/SC_results_informed$plotSigma$Mean

start_pos<-regexpr("\\_[^\\_]*$", unlist(strsplit(SC_results_informed$names, "[.]"))[c(TRUE, FALSE)])
SC_results_informed$plotN$Year<-substring(SC_results_informed$names,28,31)
SC_results_informed$plotN$`DataType`<-substring(SC_results_informed$names,start_pos+1,nchar(SC_results_informed$names)-6)
SC_results_informed$plotN$Time<-as.factor(str_split_fixed(SC_results_informed$names, "_", 6)[,5])
SC_results_informed$plotN
SC_results_informed$plotN$Time<-factor(SC_results_informed$plotN$Time,levels=c("3mos"))
SC_results_informed$plotN$Model<-rep("SC",nrow(SC_results_informed$plotN))
SC_results_informed$plotN$Prior<-rep("Gamma",nrow(SC_results_informed$plotN))
time.labs <- c( "Aug 1 - Oct 31")
names(time.labs) <- levels(SC_results_informed$plotN$Time)
SC_results_informed$plotN

SC_results_informed$plotLam0$Year<-substring(SC_results_informed$names,28,31)
SC_results_informed$plotLam0$`DataType`<-substring(SC_results_informed$names,start_pos+1,nchar(SC_results_informed$names)-6)
SC_results_informed$plotLam0$Time<-as.factor(str_split_fixed(SC_results_informed$names, "_", 6)[,5])
head(SC_results_informed$plotLam0)
SC_results_informed$plotLam0$Time<-factor(SC_results_informed$plotLam0$Time,levels=c("3mos"))
SC_results_informed$plotLam0$Model<-rep("SC",nrow(SC_results_informed$plotLam0))
SC_results_informed$plotLam0$Prior<-rep("Gamma",nrow(SC_results_informed$plotLam0))
time.labs <- c( "Aug 1 - Oct 31")
names(time.labs) <- levels(SC_results_informed$plotLam0$Time)
SC_results_informed$plotLam0

SC_results_informed$plotSigma$Year<-substring(SC_results_informed$names,28,31)
SC_results_informed$plotSigma$`DataType`<-substring(SC_results_informed$names,start_pos+1,nchar(SC_results_informed$names)-6)
SC_results_informed$plotSigma$Time<-as.factor(str_split_fixed(SC_results_informed$names, "_", 6)[,5])
head(SC_results_informed$plotSigma)
SC_results_informed$plotSigma$Time<-factor(SC_results_informed$plotSigma$Time,levels=c("Spring","4mos","3mos","Winter"))
SC_results_informed$plotSigma$Model<-rep("SC",nrow(SC_results_informed$plotSigma))
SC_results_informed$plotSigma$Prior<-rep("Gamma",nrow(SC_results_informed$plotSigma))
time.labs <- c( "Spring: May 1 - Jul 15","Jul 1 - Oct 31", "Aug 1 - Oct 31", "Winter: Sep 15 - Dec 15")
names(time.labs) <- levels(SC_results_informed$plotSigma$Time)
SC_results_informed$plotSigma


#Informed SPIM results
load("output_SPIMresults_Informative.RData")
SPIM_results_informed<-tmp_list
SPIM_results_informed$plotN$CoV<-SPIM_results_informed$plotN$SD/SPIM_results_informed$plotN$Mean
SPIM_results_informed$plotn$CoV<-SPIM_results_informed$plotn$SD/SPIM_results_informed$plotn$Mean
SPIM_results_informed$plotLam0$CoV<-SPIM_results_informed$plotLam0$SD/SPIM_results_informed$plotLam0$Mean
SPIM_results_informed$plotSigma$CoV<-SPIM_results_informed$plotSigma$SD/SPIM_results_informed$plotSigma$Mean

names<-do.call(rbind.data.frame, str_split(SPIM_results_informed$names,"_"))
colnames(names)<-seq(1,ncol(names),1)
SPIM_results_informed$plotN$Model<-rep("SPIM",nrow(SPIM_results_informed$plotN))
SPIM_results_informed$plotN$Prior<-rep("Gamma",nrow(SPIM_results_informed$plotN))
SPIM_results_informed$plotN$DataType<-names[,5]
SPIM_results_informed$plotN$Year<-names[,4]
SPIM_results_informed$plotN

SPIM_results_informed$plotLam0$Year<-SPIM_results_informed$plotN$Year
SPIM_results_informed$plotLam0$DataType<-SPIM_results_informed$plotN$DataType
SPIM_results_informed$plotLam0$Model<-rep("SPIM",nrow(SPIM_results_informed$plotLam0))
SPIM_results_informed$plotLam0$Prior<-rep("Gamma",nrow(SPIM_results_informed$plotLam0))
SPIM_results_informed$plotLam0

SPIM_results_informed$plotSigma$Year<-SPIM_results_informed$plotN$Year
SPIM_results_informed$plotSigma$DataType<-SPIM_results_informed$plotN$DataType
SPIM_results_informed$plotSigma$Model<-rep("SPIM",nrow(SPIM_results_informed$plotSigma))
SPIM_results_informed$plotSigma$Prior<-rep("Gamma",nrow(SPIM_results_informed$plotSigma))
SPIM_results_informed$plotSigma

SPIM_results_informed$plotn$Year<-SPIM_results_informed$plotN$Year
SPIM_results_informed$plotn$DataType<-SPIM_results_informed$plotN$DataType
SPIM_results_informed$plotn$Model<-rep("SPIM",nrow(SPIM_results_informed$plotn))
SPIM_results_informed$plotn$Prior<-rep("Gamma",nrow(SPIM_results_informed$plotn))
SPIM_results_informed$plotn



###plots
### SC vs SPIM, informed vs uninformed plots #####

#plot SC informed vs uninformed
SC_resultsN<-rbind(SC_results_informed$plotN[SC_results_informed$plotN$DataType=="sum"&SC_results_informed$plotN$Time=="3mos",-which(colnames(SC_results_informed$plotN)=="Time")],
                        SC_results_uninformed$plotN[SC_results_uninformed$plotN$Time=="3mos",-which(colnames(SC_results_uninformed$plotN)=="Time")])
ggplot(SC_resultsN)+
  geom_segment(  aes(x=Prior,y=HPDI_5,xend=Prior, yend=HPDI_95,color=Prior), size=1)+
  geom_segment(  aes(x=Prior, y=HPDI_25,xend=Prior, yend=HPDI_50,color=Prior), size=2)+
  geom_point(  aes(x=Prior, y=Mode), size=3, shape=21, fill="white") +
  geom_text_repel( nudge_x=0.15,aes(x = Prior,  y = Mode,  label = round(Mode,digits = 0)))+
  coord_cartesian(ylim=c(0, 300))+ 
  facet_grid(cols=vars(Year),switch="y",scales = "free")+ #
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+#,strip.placement = "outside")+
  labs(title="Caribou Abundance in Algar with SC Models",
       #subtitle = "Informed vs. Uninformed Sigma Prior",
       x ="", y = "Abundance estimate")    


#plot SPIM informed vs uninformed
SPIM_resultsN<-SPIM_results_informed$plotN)
SPIM_resultsN$DataType[which(SPIM_resultsN$DataType=="grpSize1")]<-"GrpSize1"
ggplot(SPIM_resultsN)+
  geom_segment(  aes(x=Prior,y=HPDI_5,xend=Prior, yend=HPDI_95,color=Prior), size=1)+
  geom_segment(  aes(x=Prior, y=HPDI_25,xend=Prior, yend=HPDI_50,color=Prior), size=2)+
  geom_point(  aes(x=Prior, y=Mode), size=3, shape=21, fill="white") +
  geom_text_repel( nudge_x=0.15,aes(x = Prior,  y = Mode,  label = round(Mode,digits = 0)))+
  coord_cartesian(ylim=c(0, 800))+ 
  facet_grid(cols=vars(Year),rows=vars(DataType),switch="y",scales = "free")+ #
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+#,strip.placement = "outside")+
  labs(title="Caribou Abundance in Algar in SPIM Models",
       #subtitle = "Informed vs. Uninformed Sigma Prior",
       x ="", y = "Abundance estimate")


#plot  N, with just "sum" for count version of SC data 
SC_SPIM_resultsN<-rbind(SC_results_informed$plotN[SC_results_informed$plotN$DataType=="sum"&SC_results_informed$plotN$Time=="3mos",-which(colnames(SC_results_informed$plotN)=="Time")],
                        SC_results_uninformed$plotN[SC_results_uninformed$plotN$Time=="3mos",-which(colnames(SC_results_uninformed$plotN)=="Time")],
                        SPIM_results_informed$plotN)#,
                        #SPIM_results_Uninformed$plotN)
SC_SPIM_resultsN$DataType[which(SC_SPIM_resultsN$DataType=="grpSize1")]<-"GrpSize1"
SC_SPIM_resultsN$DataType[which(SC_SPIM_resultsN$DataType=="sum")]<-"SC"
SC_SPIM_resultsN$DataType[which(SC_SPIM_resultsN$DataType=="GrpSize1")]<-"SPIM\nGrp Size 1"
SC_SPIM_resultsN$DataType[which(SC_SPIM_resultsN$DataType=="all")]<-"SPIM\nAll"
SC_SPIM_resultsN$DataType<-factor(SC_SPIM_resultsN$DataType,levels=c("SC","SPIM\nGrp Size 1","SPIM\nAll"))
SC_SPIM_resultsN$Prior<-factor(SC_SPIM_resultsN$Prior,levels=c("Gamma","Uniform"))


#add Density on right 
area<-1505
SC_SPIM_resultsD<-SC_SPIM_resultsN
SC_SPIM_resultsD[,1:7]<-SC_SPIM_resultsD[,1:7]/area*1000
SC_SPIM_resultsN$parm<-"N"
SC_SPIM_resultsD$parm<-"D"
SC_SPIM_results_ND<-rbind(SC_SPIM_resultsN,SC_SPIM_resultsD)

ggplot(SC_SPIM_resultsN)+
  geom_segment(  aes(x=DataType,y=HPDI_5,xend=DataType, yend=HPDI_95,color=Model), size=1)+
  geom_segment(  aes(x=DataType, y=HPDI_25,xend=DataType, yend=HPDI_50,color=Model), size=2)+
  geom_point(  aes(x=DataType, y=Mode), size=3, shape=21, fill="white") +
  geom_text_repel( nudge_x=0.05,aes(x = DataType,  y = Mode,  label = round(Mode,digits = 0)))+
  coord_cartesian(ylim=c(0, 800))+ 
  facet_grid(cols=vars(Year),rows=vars(Prior),switch="y",scales = "free")+ #
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+#,strip.placement = "outside")+
  labs(title="Caribou Abundance in Algar: SC vs. SPIM Models",
       #subtitle = "Informed vs. Uninformed Sigma Prior",
       x ="Data Type ", y = "Abundance estimate")+
  scale_y_continuous(sec.axis = sec_axis(~ . /area*1000, name = "Density per 1000 km2"))

#or do density in left y axis and N on right y axis
#with McFarlane et al. 2020's ESAR estimate
ggplot(SC_SPIM_resultsD, aes(x=as.numeric(DataType)))+
  geom_hline(yintercept=50.6,lty=2,color="darkgray")+
  geom_ribbon(aes(x=as.numeric(DataType), ymin=42.9, ymax=59.6),alpha=0.20)+
  geom_segment(  aes(x=as.numeric(DataType),y=HPDI_5,xend=as.numeric(DataType), yend=HPDI_95,color=Model), size=1)+
  geom_segment(  aes(x=as.numeric(DataType), y=HPDI_25,xend=as.numeric(DataType), yend=HPDI_50,color=Model), size=2)+
  geom_point(  aes(x=as.numeric(DataType), y=Mode), size=3, shape=21, fill="white") +
  geom_text_repel(min.segment.length = Inf, nudge_x=0.15,aes(x = as.numeric(DataType),  y = Mode,  label = round(Mode,digits = 1)))+
  coord_cartesian(ylim=c(0, 600))+ 
  facet_grid(Prior~Year,switch="y")+ #
  theme(panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none",plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title="Caribou Density and Abundance in Algar",#\nSC vs. SPIM Models",
       x ="Model ", y = "Density (N/1000 km2)")+
  scale_x_continuous(breaks = 1:3, labels = levels(SC_SPIM_resultsD$DataType))+
  scale_y_continuous(sec.axis = sec_axis(~ . *area/1000, name = "Abundance (N)"))

#just the spim results with all data an the gamma prior
temp<-rbind(SC_SPIM_resultsD[SC_SPIM_resultsD$Prior=="Gamma" & SC_SPIM_resultsD$DataType=="SPIM\nAll",],
SC_SPIM_resultsD[SC_SPIM_resultsD$Prior=="Gamma" & SC_SPIM_resultsD$DataType=="SC",])

nudge_amount<-c(-0.05,0.05)
nudge_direction<-c(-1,1)

plot_basic<-ggplot(temp,aes(x=as.numeric(Year)))+
  geom_segment(aes(x = 2015.5, y = 50.6, xend = 2019.6, yend = 50.6), lty=2)+
  geom_ribbon(aes(x=c(2015.5, as.numeric(Year[-c(1,2)]),2019.6), ymin=42.9, ymax=59.6),alpha=0.25)+
  geom_segment(  aes(x=as.numeric(Year),y=HPDI_5,xend=as.numeric(Year), yend=HPDI_95,color=Model), position = position_nudge(x = nudge_amount[as.numeric(as.factor(temp$Model))], y = 0), size=1)+
  geom_segment(  aes(x=as.numeric(Year), y=HPDI_25,xend=as.numeric(Year), yend=HPDI_50,color=Model),  position = position_nudge(x = nudge_amount[as.numeric(as.factor(temp$Model))], y = 0),size=2)+
  geom_point(  aes(x=as.numeric(Year), y=Mode), position = position_nudge(x = nudge_amount[as.numeric(as.factor(temp$Model))], y = 0) ,size=3, shape=21, fill="white") +
  geom_text_repel(min.segment.length = Inf,
                  nudge_x=0.2*nudge_direction[as.numeric(as.factor(temp$Model))],
                  aes(x = as.numeric(Year),  y = Mode,  label = round(Mode,digits = 1)),size=5)+
  theme(panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none",plot.title = element_text(size=20,hjust = 0.5),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=15,angle = 0, hjust = 0.5),
        axis.text.y = element_text(size=15,angle = 0, hjust = 0.5),
        axis.title = element_text(size = 12),
        plot.subtitle = element_text(hjust = 0.5),)+#,strip.placement = "outside")+
  labs(title="Algar Caribou Density and Abundance\nSC vs. SPIM Models",#, \n SPIM (All Data ; Gamma Prior)
       x ="Year ", y = "Density (N per 1000 km2)")+
  coord_cartesian(xlim = c(2015.7,2019.4))+
  scale_y_continuous(sec.axis = sec_axis(~ . *area/1000, name = "Abundance (N)"))


nudge_amount<-c(-0.05,0.05)
nudge_direction<-c(-1,1)  

plot_abundance<-ggplot(SC_SPIM_resultsN)+
  geom_segment(  aes(x=Year,y=HPDI_5,xend=Year, yend=HPDI_95,color=Study), position = position_nudge(x = nudge_amount[as.numeric(as.factor(tmp_N$Study))], y = 0), size=1)+
  geom_segment(  aes(x=Year, y=HPDI_25,xend=Year, yend=HPDI_50,color=Study), position = position_nudge(x = nudge_amount[as.numeric(as.factor(tmp_N$Study))], y = 0), size=2)+
  geom_point(  aes(x=Year, y=Mode),  position = position_nudge(x = nudge_amount[as.numeric(as.factor(tmp_N$Study))], y = 0),size=3, shape=21, fill="white") +
  #geom_text_repel(position = position_nudge(x = nudge_more[as.numeric(as.factor(tmp_N$Model))], y = 0),segment.color = NA,aes(x = Year,  y = Mode,  label = round(Mode,digits = 1)))+
  geom_text_repel(nudge_x=0.17*nudge_direction[as.numeric(as.factor(tmp_N$Study))],segment.color = NA,aes(x = Year,  y = Mode,  label = round(Mode,digits = 1)))+
  coord_cartesian(ylim=c(0, 800))+ 
  facet_grid(cols=vars(Model),switch="y",scales = "free")+ #
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+#,strip.placement = "outside")+
  labs(title="Caribou Abundance across Landscapes",
       x ="", y = "Abundance estimate") #Density ( N/100 km2)
plot_abundance


#plot sigma
SC_SPIM_resultsSigma<-rbind(SC_results_informed$plotSigma[SC_results_informed$plotSigma$DataType=="sum"&SC_results_informed$plotSigma$Time=="3mos",-which(colnames(SC_results_informed$plotSigma)=="Time")],
                        SC_results_uninformed$plotSigma[SC_results_uninformed$plotSigma$Time=="3mos",-which(colnames(SC_results_uninformed$plotSigma)=="Time")],
                        SPIM_results_informed$plotSigma)#,
                        #SPIM_results_Uninformed$plotSigma)
SC_SPIM_resultsSigma$DataType[which(SC_SPIM_resultsSigma$DataType=="grpSize1")]<-"GrpSize1"
SC_SPIM_resultsSigma$DataType[which(SC_SPIM_resultsSigma$DataType=="sum")]<-"SC"
SC_SPIM_resultsSigma$DataType[which(SC_SPIM_resultsSigma$DataType=="GrpSize1")]<-"SPIM\nGrp Size 1"
SC_SPIM_resultsSigma$DataType[which(SC_SPIM_resultsSigma$DataType=="all")]<-"SPIM\nAll"
SC_SPIM_resultsSigma$DataType<-factor(SC_SPIM_resultsSigma$DataType,levels=c("SC","SPIM\nGrp Size 1","SPIM\nAll"))
SC_SPIM_resultsSigma$Prior<-factor(SC_SPIM_resultsSigma$Prior,levels=c("Gamma","Uniform"))

ggplot(SC_SPIM_resultsSigma)+
  geom_segment(  aes(x=DataType,y=HPDI_5,xend=DataType, yend=HPDI_95,color=Model), size=1)+
  geom_segment(  aes(x=DataType, y=HPDI_25,xend=DataType, yend=HPDI_50,color=Model), size=2)+
  geom_point(  aes(x=DataType, y=Mode), size=3, shape=21, fill="white") +
  geom_text_repel(min.segment.length = Inf, nudge_x=0.15,aes(x = DataType,  y = Mode,  label = round(Mode,digits = 0)))+
  coord_cartesian(ylim=c(0, 10))+ 
  facet_grid(cols=vars(Year),rows=vars(Prior),switch="y",scales = "free")+ #
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+#,strip.placement = "outside")+
  labs(title="Caribou Sigma in Algar: SC vs. SPIM Models",
       x ="Data Type ", y = "Sigma estimate (km)")    


ggplot(data=SPIM_results$plotn)+
  geom_segment(  aes(x=DataType,y=HPDI_5,xend=DataType, yend=HPDI_95), size=1, color="lightblue")+
  geom_segment(  aes(x=DataType, y=HPDI_25,xend=DataType, yend=HPDI_50), size=2, color="darkblue")+
  geom_point( aes(x=DataType, y=Mode), size=3, shape=21, fill="white") +
  geom_text_repel(nudge_x=0.15,aes(x = DataType,  y = Mode,  label = round(Mode,digits = 0)))+
  coord_cartesian(ylim=c(0, 50))+ 
  facet_grid(cols=vars(Year),scales = "free_x")+ #
  theme(plot.title = element_text(hjust = 0.5))+#,strip.placement = "outside")+
  labs(title="Estimated number of observed individuals using Spatial Partial ID Models",
       x ="Data type", y = "Observed Individuals")
