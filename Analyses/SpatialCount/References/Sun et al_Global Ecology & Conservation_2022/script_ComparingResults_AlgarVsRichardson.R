#combining/comparing SC and SPIm results for caribou
#across landscapes!
library(ggplot2)
library(stringr)
library(ggrepel)
library(ggpubr)
library(rphylopic)
library(grid)

#####Go Run the script for Richardson, and then 
#sc informed
load("output_SCresults_informativePriors.RData")
SC_results_informed_Rich<-tmp_list
SC_results_informed_Rich_names<-str_split_fixed(SC_results_informed_Rich$names,"_",7)
SC_results_informed_Rich$plotN$CoV<-SC_results_informed_Rich$plotN$SD/SC_results_informed_Rich$plotN$Mean
SC_results_informed_Rich$plotSigma$CoV<-SC_results_informed_Rich$plotSigma$SD/SC_results_informed_Rich$plotSigma$Mean
SC_results_informed_Rich$plotLam0$CoV<-SC_results_informed_Rich$plotLam0$SD/SC_results_informed_Rich$plotLam0$Mean

#spim informed
load("output_SPIMresults_Informative.RData")
SPIM_results_informed_Rich<-tmp_list
SPIM_results_informed_Rich_names<-str_split_fixed(SPIM_results_informed_Rich$names, "_", 7)
SPIM_results_informed_Rich$plotN$CoV<-SPIM_results_informed_Rich$plotN$SD/SPIM_results_informed_Rich$plotN$Mean
SPIM_results_informed_Rich$plotn$CoV<-SPIM_results_informed_Rich$plotn$SD/SPIM_results_informed_Rich$plotn$Mean
SPIM_results_informed_Rich$plotLam0$CoV<-SPIM_results_informed_Rich$plotLam0$SD/SPIM_results_informed_Rich$plotLam0$Mean
SPIM_results_informed_Rich$plotSigma$CoV<-SPIM_results_informed_Rich$plotSigma$SD/SPIM_results_informed_Rich$plotSigma$Mean
SPIM_results_informed_Rich

#rename to add "Rich", Year, Time, Prior, DataType
#this starts with 2 because 1 is just the vector of model output names
for(i in 2:4){ #to plotN, plotLam0, and plotSigma
  SC_results_informed_Rich[[i]]$Study<-"Richardson"
  SPIM_results_informed_Rich[[i]]$Study<-"Richardson"
  SPIM_results_informed_Rich[[i+1]]$Study<-"Richardson"
  
  SC_results_informed_Rich[[i]]$Prior<-"Gamma"
  SPIM_results_informed_Rich[[i]]$Prior<-"Gamma"
  SPIM_results_informed_Rich[[i+1]]$Prior<-"Gamma"
  
  SC_results_informed_Rich[[i]]$Time<-SC_results_informed_Rich_names[,5]
  SPIM_results_informed_Rich[[i]]$Time<-"3mos"
  SPIM_results_informed_Rich[[i+1]]$Time<-"3mos"
  
  SC_results_informed_Rich[[i]]$Year<-SC_results_informed_Rich_names[,4]
  SPIM_results_informed_Rich[[i]]$Year<-SPIM_results_informed_Rich_names[,4]
  SPIM_results_informed_Rich[[i+1]]$Year<-SPIM_results_informed_Rich_names[,4]
  
  SC_results_informed_Rich[[i]]$DataType<-substring(SC_results_informed_Rich_names[,6],1,str_locate(SC_results_informed_Rich_names[,6], "\\.")[,1]-1)
  SPIM_results_informed_Rich[[i]]$DataType<-SPIM_results_informed_Rich_names[,5]
  SPIM_results_informed_Rich[[i+1]]$DataType<-SPIM_results_informed_Rich_names[,5]
  SPIM_results_informed_Rich[[i]]$DataType[which(SPIM_results_informed_Rich_names[,6]=="wCollars")]<- "allWCollars"                                           
  SPIM_results_informed_Rich[[i+1]]$DataType[which(SPIM_results_informed_Rich_names[,6]=="wCollars")]<-"allWCollars"
}

SC_results_informed_Rich$plotN
SPIM_results_informed_Rich$plotN


####Go Run Run the Scripts for Algar and then
load("output_SCresults_informativePriors.RData")
SC_results_informed_Algar<-tmp_list
SC_results_informed_Algar_names<-str_split_fixed(SC_results_informed_Algar$names, "_", 7)
SC_results_informed_Algar$plotN$CoV<-SC_results_informed_Algar$plotN$SD/SC_results_informed_Algar$plotN$Mean
SC_results_informed_Algar$plotSigma$CoV<-SC_results_informed_Algar$plotSigma$SD/SC_results_informed_Algar$plotSigma$Mean
SC_results_informed_Algar$plotLam0$CoV<-SC_results_informed_Algar$plotLam0$SD/SC_results_informed_Algar$plotLam0$Mean
SC_results_informed_Algar

load("output_SPIMresults_Informative.RData")
SPIM_results_informed_Algar<-tmp_list
SPIM_results_informed_Algar_names<-str_split_fixed(SPIM_results_informed_Algar$names, "_", 7)
SPIM_results_informed_Algar$plotN$CoV<-SPIM_results_informed_Algar$plotN$SD/SPIM_results_informed_Algar$plotN$Mean
SPIM_results_informed_Algar$plotn$CoV<-SPIM_results_informed_Algar$plotn$SD/SPIM_results_informed_Algar$plotn$Mean
SPIM_results_informed_Algar$plotLam0$CoV<-SPIM_results_informed_Algar$plotLam0$SD/SPIM_results_informed_Algar$plotLam0$Mean
SPIM_results_informed_Algar$plotSigma$CoV<-SPIM_results_informed_Algar$plotSigma$SD/SPIM_results_informed_Algar$plotSigma$Mean
SPIM_results_informed_Algar

#rename to add "Algar", Year, Time, Prior, DataType
for(i in 2:4){ #for plotN, plotLam0, and plotSigma
  SC_results_informed_Algar[[i]]$Study<-"Algar"
  SPIM_results_informed_Algar[[i]]$Study<-"Algar"
  SPIM_results_informed_Algar[[i+1]]$Study<-"Algar"
  
  SC_results_informed_Algar[[i]]$Time<-SC_results_informed_Algar_names[,5]
  SPIM_results_informed_Algar[[i]]$Time<-"3mos"
  SPIM_results_informed_Algar[[i+1]]$Time<-"3mos"
  
  SC_results_informed_Algar[[i]]$Prior<-"Gamma"
  SPIM_results_informed_Algar[[i]]$Prior<-"Gamma"
  SPIM_results_informed_Algar[[i+1]]$Prior<-"Gamma"
  
  SC_results_informed_Algar[[i]]$Year<-SC_results_informed_Algar_names[,4]
  SPIM_results_informed_Algar[[i]]$Year<-SPIM_results_informed_Algar_names[,4]
  SPIM_results_informed_Algar[[i+1]]$Year<-SPIM_results_informed_Algar_names[,4]
  

  SC_results_informed_Algar[[i]]$DataType<-substring(SC_results_informed_Algar_names[,6],1,str_locate(SC_results_informed_Algar_names[,6], "\\.")[,1]-1)
  
  SPIM_results_informed_Algar[[i]]$DataType<-as.character(SPIM_results_informed_Algar[[i]]$DataType)
  SPIM_results_informed_Algar[[i]]$DataType[which(SPIM_results_Uninformed_Algar_names[,6]=="wCollars")]<- "allWCollars"                                           
  SPIM_results_informed_Algar[[i+1]]$DataType<-as.character(SPIM_results_informed_Algar[[i+1]]$DataType)
  SPIM_results_informed_Algar[[i+1]]$DataType[which(SPIM_results_informed_Algar_names[,6]=="wCollars")]<- "allWCollars"                                           
                                            
  
}

SC_results_informed_Algar$plotN
SPIM_results_informed_Algar$plotN


###plot all the N data, as is
all_N<-rbind(SC_results_informed_Algar$plotN,
             SC_results_informed_Rich$plotN,
             SPIM_results_informed_Algar$plotN,
             SPIM_results_informed_Rich$plotN)

all_N$Model<-c(rep("SC",nrow(SC_results_informed_Algar$plotN)+
                     nrow(SC_results_informed_Rich$plotN)),
               rep("SPIM",nrow(SPIM_results_informed_Algar$plotN)+
                     nrow(SPIM_results_informed_Rich$plotN)))

all_N$DataType[all_N$DataType%in%c("allWCollars", "allwCollar" )]<-"all"
all_N$DataType<-factor(all_N$DataType,levels = c("sum", "grpSize1", "all"))

#####

###plot all the sigma data, as is
all_Sigma<-rbind(SC_results_informed_Algar$plotSigma,
             SC_results_informed_Rich$plotSigma,
             SPIM_results_informed_Algar$plotSigma,
             SPIM_results_informed_Rich$plotSigma)

all_Sigma$Model<-c(rep("SC",nrow(SC_results_informed_Algar$plotSigma)+
                     nrow(SC_results_informed_Rich$plotSigma)),
               rep("SPIM",nrow(SPIM_results_informed_Algar$plotSigma)+
                     nrow(SPIM_results_informed_Rich$plotSigma)))
all_Sigma$DataType[all_N$DataType%in%c("allWCollars", "allwCollar" )]<-"all"
all_Sigma$DataType<-factor(all_N$DataType,levels = c("sum", "grpSize1", "all"))


###Now for some more fine tuning and parsing


#combine informed SC results across different landscapes
SC_ResultsN_axLand<-rbind(SC_results_informed_Rich$plotN[SC_results_informed_Rich$plotN$Time=="3mos",],
                          SC_results_informed_Algar$plotN[SC_results_informed_Algar$plotN$Time=="3mos",])

SC_ResultsSigma_axLand<-rbind(SC_results_informed_Rich$plotSigma[SC_results_informed_Rich$plotSigma$Time=="3mos",],
                              SC_results_informed_Algar$plotSigma[SC_results_informed_Algar$plotSigma$Time=="3mos",])

SC_ResultsLam_axLand<-rbind(SC_results_informed_Rich$plotLam0[SC_results_informed_Rich$plotLam0$Time=="3mos",],
                            SC_results_informed_Algar$plotLam0[SC_results_informed_Algar$plotLam0$Time=="3mos",] )

#combine informedSPIM results across different landscapes
SPIM_ResultsN_axLand<-rbind(SPIM_results_informed_Rich$plotN,
                            SPIM_results_informed_Algar$plotN)
 
SPIM_ResultsSigma_axLand<-rbind(SPIM_results_informed_Rich$plotSigma,
                               SPIM_results_informed_Algar$plotSigma)
 
SPIM_ResultsLam_axLand<-rbind(SPIM_results_informed_Rich$plotLam0,
                              SPIM_results_informed_Algar$plotLam0)

#combine informedSPIM and SC results across different landscapes 
SC_SPIM_ResultsN_axLand<-rbind(SC_ResultsN_axLand,SPIM_ResultsN_axLand)
SC_SPIM_ResultsN_axLand$Model<-c(rep("SC",nrow(SC_ResultsN_axLand)),rep("SPIM",nrow(SPIM_ResultsN_axLand)))



########################
####  D  plots #########
########################

#combine informed SPIM and SC results across different landscapes 
#for some of the below plots, will need to reinstate the !=GrpSize1
SC_SPIM_ResultsN_axLand<-rbind(SC_ResultsN_axLand[SC_ResultsN_axLand$DataType=="sum",],
                               SPIM_ResultsN_axLand)
SC_SPIM_ResultsN_axLand$Model<-c(rep("SC",nrow(SC_ResultsN_axLand[SC_ResultsN_axLand$DataType=="sum",])),#
                                 rep("SPIM",nrow(SPIM_ResultsN_axLand)))

#plot informed SC and SPIM  N's
nudge_amount<-c(-0.05,0.05)
nudge_direction<-c(-1,1)
as.numeric(as.factor(SC_SPIM_ResultsN_axLand$Study))

plot_abundance<-ggplot(SC_SPIM_ResultsN_axLand)+ #data=tmp_N
  geom_segment(  aes(x=Year,y=HPDI_5,xend=Year, yend=HPDI_95,color=Study), position = position_nudge(x = nudge_amount[as.numeric(as.factor(SC_SPIM_ResultsN_axLand$Study))], y = 0), size=1)+
  geom_segment(  aes(x=Year, y=HPDI_25,xend=Year, yend=HPDI_50,color=Study), position = position_nudge(x = nudge_amount[as.numeric(as.factor(SC_SPIM_ResultsN_axLand$Study))], y = 0), size=2)+
  geom_point(  aes(x=Year, y=Mode),  position = position_nudge(x = nudge_amount[as.numeric(as.factor(SC_SPIM_ResultsN_axLand$Study))], y = 0),size=3, shape=21, fill="white") +
  #geom_text_repel(position = position_nudge(x = nudge_more[as.numeric(as.factor(SC_SPIM_ResultsN_axLand$Model))], y = 0),segment.color = NA,aes(x = Year,  y = Mode,  label = round(Mode,digits = 1)))+
  geom_text_repel(nudge_x=0.17*nudge_direction[as.numeric(as.factor(SC_SPIM_ResultsN_axLand$Study))],segment.color = NA,aes(x = Year,  y = Mode,  label = round(Mode,digits = 1)))+
  coord_cartesian(ylim=c(0, 800))+ 
  facet_grid(cols=vars(Model),switch="y",scales = "free")+ #
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  labs(title="Caribou Abundance across Landscapes\n(all data w collars)",
       x ="", y = "Abundance estimate") #Density ( N/100 km2)
plot_abundance

#plot DENSITY, with informed priors, just "sum" for count version of SC data 
AlgarArea<-1504.934
RichArea<-1999

SC_SPIM_ResultsD_axLand<- SC_SPIM_ResultsN_axLand
SC_SPIM_ResultsD_axLand$Nmode<-SC_SPIM_ResultsD_axLand$Mode
SC_SPIM_ResultsD_axLand[which(SC_SPIM_ResultsD_axLand$Study=="Algar"),1:8]<-SC_SPIM_ResultsD_axLand[which(SC_SPIM_ResultsD_axLand$Study=="Algar"),1:8]/AlgarArea*1000
SC_SPIM_ResultsD_axLand[which(SC_SPIM_ResultsD_axLand$Study=="Richardson"),1:8]<-SC_SPIM_ResultsD_axLand[which(SC_SPIM_ResultsD_axLand$Study=="Richardson"),1:8]/RichArea*1000


#### density visualization ####

scSPIM_Algar<-SC_SPIM_ResultsD_axLand[SC_SPIM_ResultsD_axLand$Study=="Algar"&
                                        SC_SPIM_ResultsD_axLand$DataType!="grpSize1",]
scSPIM_Algar$Year<-as.factor(scSPIM_Algar$Year)
nudge_order<-as.numeric(as.factor(scSPIM_Algar$Model))
nudge_amount3<-c(0,0.1)
color_order<-c("blue","red")
shape_order<-c(21,16)

scSPIM_Algar$color<-color_order[nudge_order]
scSPIM_Algar$`segment legend`<-"McFarlane et al. 2020"
scSPIM_Algar$`segment legend2`<-"Minimum"
plot_Algar_density<-ggplot(data=scSPIM_Algar,
                        aes(x= as.numeric(as.character(Year))))+ 
   geom_segment(aes(x = 2015.5, y = 50.6, xend = 2019.6, yend = 50.6,linetype=`segment legend`),
                color="darkgray",show.legend = FALSE)+
  geom_segment(aes(x = 2015.5, y = 17.2, xend = 2019.6, yend = 17.2,linetype=`segment legend2`),
               lwd=1,color="darkgray",show.legend = FALSE)+
  geom_ribbon(aes(x=c(2015.5,2019.56,2019.57,2019.58, 2019.59 ,2019.598 ,2019.599 ,2019.6),
                   ymin=42.9, ymax=59.6),alpha=0.20)+
  geom_linerange(aes( ymin=HPDI_5,ymax=HPDI_95, color=Model,group=Model),
                 position=position_nudge(nudge_amount3[nudge_order]),size=1,show.legend = FALSE)+
  geom_linerange( aes( ymin=HPDI_25, ymax=HPDI_50,color=Model,group=Model),
                  position=position_nudge(nudge_amount3[nudge_order]),size=2,show.legend = FALSE)+
  geom_point( aes( y=Mode,shape=Model),
              position=position_nudge(nudge_amount3[nudge_order]),size=2, 
              color="black",fill="white",show.legend = FALSE)+
  coord_cartesian(ylim=c(0, 400), xlim=c(2015.7,2019.4))+
  labs(title="Algar",
       x ="",y="")+# y =expression(paste("Density (N / 1000 ", km^2,")")))+
  theme(axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_text(size=20,hjust = 0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA))+
  scale_shape_manual("Model", values=shape_order)+
  scale_color_manual(name="Model",
                     values=c("blue","red"),
                     labels=c("SC","SPIM"))+
  scale_linetype_manual("",values=c("McFarlane et al. 2020"=2,
                                    "Minimum"=5))


scSPIM_Rich<-SC_SPIM_ResultsD_axLand[SC_SPIM_ResultsD_axLand$Study=="Richardson"&
                                        SC_SPIM_ResultsD_axLand$DataType!="grpSize1",]
scSPIM_Rich$Year<-as.factor(scSPIM_Rich$Year)
nudge_order<-as.numeric(as.factor(scSPIM_Rich$Model))
nudge_amount3<-c(0,0.1)
color_order<-c("blue","red")
shape_order<-c(21,16)

scSPIM_Rich$color<-color_order[nudge_order]
scSPIM_Rich$`segment legend`<-"Minimum"
plot_Rich_density<-ggplot(data=scSPIM_Rich,
                           aes(x= as.numeric(as.character(Year))))+ #data=tmp_N
  geom_segment(aes(x = 2015.5, y = 17.7, xend = 2019.6, yend = 17.7,linetype=`segment legend`),
               lwd=1,color="darkgray",show.legend = FALSE)+
  geom_linerange(aes( ymin=HPDI_5,ymax=HPDI_95, color=Model,group=Model),
                 position=position_nudge(nudge_amount3[nudge_order]),size=1,show.legend = FALSE)+
  geom_linerange( aes( ymin=HPDI_25, ymax=HPDI_50,color=Model,group=Model),
                  position=position_nudge(nudge_amount3[nudge_order]),size=2,show.legend = FALSE)+
  geom_point( aes( y=Mode,shape=Model),
              position=position_nudge(nudge_amount3[nudge_order]),size=2, 
              color="black",fill="white",show.legend = FALSE)+
  coord_cartesian(ylim=c(0, 400), xlim=c(2015.7,2019.4))+
  labs(title="Richardson",
       x ="", y =expression(paste("Density (N / 1000 ", km^2,")")))+
  theme( axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_text(size=20,hjust = 0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA))+
  scale_shape_manual("Model", values=shape_order)+
  scale_color_manual(name="Model",
                     values=c("blue","red"),
                     labels=c("SC","SPIM"))+
  scale_linetype_manual("",values=c("Minimum"=5))


plot_density_AlgarRichv2<-ggpubr::ggarrange(plot_Rich_density,plot_Algar_density,
                                          common.legend = TRUE, 
                                          nrow=1)
plot_density_AlgarRichv2

#############################################################################
###### density plots in the same style as plot_density_AlgarRichv2, but for ####
###### group size comparison###


tmp<-spim_AlgarRich[spim_AlgarRich$Study=="Algar",]
nudge_order99<-as.numeric((as.factor(tmp$DataType)))
nudge_amount99<-c(0,0.15)
plot_SP_Algardensity<- ggplot(data=tmp,
                        aes(x= as.numeric(as.character(Year))))+ #data=tmp_N
  
  geom_linerange(aes(ymin=HPDI_5,ymax=HPDI_95,alpha=DataType),#,alpha=alpha),
                 color="blue",
                 position=position_nudge(nudge_amount99[nudge_order99]),
                 size=1,show.legend = FALSE)+
  
  geom_linerange(aes( ymin=HPDI_25, ymax=HPDI_50,alpha=DataType), 
                 color="blue",
                 position=position_nudge(nudge_amount99[nudge_order99]),
                 size=2,show.legend = FALSE)+
  
  geom_point(aes(y=Mode),size=2,color="black",
             shape=21,
             position=position_nudge(nudge_amount99[nudge_order99]),
             fill="white",show.legend = FALSE)+
  
  labs(title="Algar",
    x ="",
    y="")+#, y = "Density (N/1000km2)"
  theme(axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_text(size=20, hjust = 0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA))+
  scale_alpha_manual(name="Detections Used",
                     values = c( 1, 0.4),labels=c("All", "Group Size = 1"))

tmp<-spim_AlgarRich[spim_AlgarRich$Study=="Richardson",]
nudge_order99<-as.numeric((as.factor(tmp$DataType)))
nudge_amount99<-c(0,0.15)

plot_SP_Richdensity<- ggplot(data=tmp,
                           aes(x= as.numeric(as.character(Year))))+ #data=tmp_N
  
  geom_linerange(aes(ymin=HPDI_5,ymax=HPDI_95,alpha=DataType),#,alpha=alpha),
                 color="red",
                 position=position_nudge(nudge_amount99[nudge_order99]),
                 size=1,show.legend = FALSE)+
  
  geom_linerange(aes( ymin=HPDI_25, ymax=HPDI_50,alpha=DataType), 
                 color="red",
                 position=position_nudge(nudge_amount99[nudge_order99]),
                 size=2,show.legend = FALSE)+
  
  geom_point(aes(y=Mode),size=2,color="black",
             shape=21,
             position=position_nudge(nudge_amount99[nudge_order99]),
             fill="white",show.legend = FALSE)+
  
  coord_cartesian(ylim=c(0, 400), xlim=c(2015.7,2019.4))+
  labs(title="Richardson",
    x ="", 
    y = expression(paste("Density (N / 1000 ", km^2,")")))+
  theme(axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_text(size=20, hjust = 0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA))+
  scale_alpha_manual(name="Detections Used",
                     values = c( 1, 0.4),labels=c("All", "Group Size = 1"))

plot_SP_densityv2<-ggpubr::ggarrange(plot_SP_Richdensity,plot_SP_Algardensity,
                                                    common.legend = TRUE,legend = "bottom", 
                                                    nrow=1)


########################
####sigma plots#########
########################
#plot sigma, with informed priors, just "sum" for count version of SC data 
SC_SPIM_ResultsSigma_axLand<-rbind(SC_ResultsSigma_axLand[SC_ResultsSigma_axLand$DataType=="sum",],
                                   SPIM_ResultsSigma_axLand)
SC_SPIM_ResultsSigma_axLand$Model<-c(rep("SC",nrow(SC_ResultsSigma_axLand[SC_ResultsSigma_axLand$DataType=="sum",])),
                                 rep("SPIM",nrow(SPIM_ResultsSigma_axLand)))


### sigma visualization####
scSPIM_Algar_sig<-SC_SPIM_ResultsSigma_axLand[SC_SPIM_ResultsSigma_axLand$Study=="Algar"&
                                                SC_SPIM_ResultsSigma_axLand$DataType!="grpSize1",]
scSPIM_Algar_sig$Year<-as.factor(scSPIM_Algar_sig$Year)
nudge_order<-as.numeric(as.factor(scSPIM_Algar_sig$Model))
nudge_amount3<-c(0,0.1)
color_order<-c("blue","red")
shape_order<-c(21,16)

scSPIM_Algar_sig$color<-color_order[nudge_order]
scSPIM_Algar_sig$`segment legend`<-"McFarlane et al. 2020"

plot_Algar_sigma<-ggplot(data=scSPIM_Algar_sig,
                      aes( x= as.numeric(as.character(Year))))+ #data=tmp_N
  geom_segment(aes(x = 2015.5, y = 1.8, xend = 2019.6, yend = 1.8,linetype=`segment legend`), color="darkgray")+
  geom_ribbon(aes(x=c(2015.5,2015.6,2015.7,2019.58 ,2019.59 ,2019.598 ,2019.599 ,2019.6),
                  ymin=1.45, ymax=2.18),alpha=0.20)+
  geom_linerange(aes(ymin=HPDI_5,ymax=HPDI_95,color=Model,group=Model,),
                 position=position_nudge(nudge_amount3[nudge_order]),size=1)+
  geom_linerange( aes(ymin=HPDI_25, ymax=HPDI_50,color=Model,group=Model,),
                  position=position_nudge(nudge_amount3[nudge_order]),size=2)+
  geom_point( aes(y=Mode,shape=Model),
              position=position_nudge(nudge_amount3[nudge_order]),size=2, 
              color="black",fill="white")+
  coord_cartesian(ylim=c(0, 5), xlim=c(2015.7,2019.4))+
  labs( x ="Year", y="")+# y = "Sigma (km)")+
  theme( axis.title = element_text(size = 15),
         axis.text = element_text(size = 12),
         legend.position="bottom",
         legend.justification = "center", #or left if you plan to pair this plot with teh SPIM plot that includes grpsize1
         legend.key=element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA))+
  scale_shape_manual("Model", values=shape_order)+
  scale_color_manual(name="Model",
                     values=c("blue","red"),
                     labels=c("SC","SPIM"))+
  scale_linetype_manual("",values=c("McFarlane et al. 2020"=2))


scSPIM_Rich_sig<-SC_SPIM_ResultsSigma_axLand[SC_SPIM_ResultsSigma_axLand$Study=="Richardson"&
                                                SC_SPIM_ResultsSigma_axLand$DataType!="grpSize1",]
scSPIM_Rich_sig$Year<-as.factor(scSPIM_Rich_sig$Year)
nudge_order<-as.numeric(as.factor(scSPIM_Rich_sig$Model))
nudge_amount3<-c(0,0.1)
color_order<-c("blue","red")
shape_order<-c(21,16)

scSPIM_Rich_sig$color<-color_order[nudge_order]

plot_Rich_sigma<-ggplot(data=scSPIM_Rich_sig,
                         aes( x= as.numeric(as.character(Year))))+ #data=tmp_N
  geom_linerange(aes(ymin=HPDI_5,ymax=HPDI_95,color=Model,group=Model,),
                 position=position_nudge(nudge_amount3[nudge_order]),size=1)+
  geom_linerange( aes(ymin=HPDI_25, ymax=HPDI_50,color=Model,group=Model,),
                  position=position_nudge(nudge_amount3[nudge_order]),size=2)+
  geom_point( aes(y=Mode,shape=Model),
              position=position_nudge(nudge_amount3[nudge_order]),size=2, 
              color="black",fill="white")+
  coord_cartesian(ylim=c(0, 5), xlim=c(2015.7,2019.4))+
  labs( x ="Year", y ="Sigma (\u03c3 ; km)")+ #"Sigma (km)")
  theme( axis.title = element_text(size = 15),
         axis.text = element_text(size = 12),
         legend.position="bottom",
         legend.justification = "center", #or left if you plan to pair this plot with teh SPIM plot that includes grpsize1
         legend.key=element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA))+
  scale_shape_manual("Model", values=shape_order)+
  scale_color_manual(name="Model",
                     values=c("blue","red"),
                     labels=c("SC","SPIM"))


#############################################################################
###### sigma plots in the same style as plot_density_AlgarRichv2, but for ####
###### group size comparison###


tmp<-spim_AlgarRich_sig[spim_AlgarRich_sig$Study=="Algar",]
nudge_order99<-as.numeric((as.factor(tmp$DataType)))
nudge_amount99<-c(0,0.15)
plot_SP_Algarsigma<- ggplot(data=tmp,
                              aes(x= as.numeric(as.character(Year))))+ #data=tmp_N
  
  geom_linerange(aes(ymin=HPDI_5,ymax=HPDI_95,alpha=DataType),#,alpha=alpha),
                 color="blue",
                 position=position_nudge(nudge_amount99[nudge_order99]),
                 size=1,show.legend = FALSE)+
  
  geom_linerange(aes( ymin=HPDI_25, ymax=HPDI_50,alpha=DataType), 
                 color="blue",
                 position=position_nudge(nudge_amount99[nudge_order99]),
                 size=2,show.legend = FALSE)+
  
  geom_point(aes(y=Mode),size=2,color="black",
             shape=21,
             position=position_nudge(nudge_amount99[nudge_order99]),
             fill="white",show.legend = FALSE)+
  
  coord_cartesian(ylim=c(0, 8), xlim=c(2015.7,2019.4))+
  labs( x ="Year",y="")+
  theme(axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_text(size=20, hjust = 0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA))+
  scale_alpha_manual(name="Detections Used",
                     values = c( 1, 0.4),labels=c("All", "Group Size = 1"))

tmp<-spim_AlgarRich_sig[spim_AlgarRich_sig$Study=="Richardson",]
nudge_order99<-as.numeric((as.factor(tmp$DataType)))
nudge_amount99<-c(0,0.15)

plot_SP_Richsigma<- ggplot(data=tmp,
                             aes(x= as.numeric(as.character(Year))))+ #data=tmp_N
  
  geom_linerange(aes(ymin=HPDI_5,ymax=HPDI_95,alpha=DataType),#,alpha=alpha),
                 color="red",
                 position=position_nudge(nudge_amount99[nudge_order99]),
                 size=1,show.legend = FALSE)+
  
  geom_linerange(aes( ymin=HPDI_25, ymax=HPDI_50,alpha=DataType), 
                 color="red",
                 position=position_nudge(nudge_amount99[nudge_order99]),
                 size=2,show.legend = FALSE)+
  
  geom_point(aes(y=Mode),size=2,color="black",
             shape=21,
             position=position_nudge(nudge_amount99[nudge_order99]),
             fill="white",show.legend = FALSE)+
  
  coord_cartesian(ylim=c(0,8), xlim=c(2015.7,2019.4))+
  labs(x ="Year",y = "Sigma (\u03c3 ; km)")+
  theme(axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_text(size=20, hjust = 0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA))+
  scale_alpha_manual(name="Detections Used",
                     values = c( 1, 0.4),labels=c("All", "Group Size = 1"))

plot_SP_sigmav2<-ggpubr::ggarrange(plot_SP_Richsigma,plot_SP_Algarsigma,
                                     common.legend = TRUE,legend = "bottom", 
                                     nrow=1)

########################
####bringing D and sigma plots together #########
########################
plot_sigma_AlgarRichv2<-ggpubr::ggarrange(plot_Rich_sigma,plot_Algar_sigma,
                                        common.legend = TRUE, 
                                        nrow=1)
plot_sigma_AlgarRichv2

plot_AlgarRichv2<-ggpubr::ggarrange(plot_Rich_density,plot_Algar_density,#plot_SP_density
                                  plot_Rich_sigma,plot_Algar_sigma,#plot_SP_sigma
                                  common.legend = TRUE,legend = "bottom", 
                                  nrow=2,ncol=2)

plot_AlgarRichv2


plot_AlgarRich_allVsgroupSize1v2<-ggpubr::ggarrange(plot_SP_densityv2,plot_SP_sigmav2,
                                                    common.legend = TRUE,legend = "bottom", 
                                                    nrow=2)




########################
#### COV plots #########
########################

SC_SPIM_ResultsD_axLand
ggplot(SC_SPIM_ResultsD_axLand,)+
  geom_boxplot(aes(x=Model,y=CoV,fill=Study))

SC_SPIM_ResultsD_axLand$PairNumber<-interaction(SC_SPIM_ResultsD_axLand$Year,SC_SPIM_ResultsD_axLand$Study)
SC_SPIM_ResultsSigma_axLand$PairNumber<-interaction(SC_SPIM_ResultsSigma_axLand$Year,SC_SPIM_ResultsSigma_axLand$Study)
SC_SPIM_ResultsD_axLand$Study<-factor(SC_SPIM_ResultsD_axLand$Study,
                                      levels=c("Richardson","Algar"))

plot_D_CoV<-ggplot(SC_SPIM_ResultsD_axLand[SC_SPIM_ResultsD_axLand$DataType!="grpSize1",], 
       aes(x = Model, y = CoV*100,label=Year)) +
  geom_line(aes(group = PairNumber, color=Year),lwd=1.5,alpha=0.5)+
  geom_point(size=4,aes(fill=Year),pch=21, color="black")+
 # geom_text()+
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 20))+
  facet_wrap(vars(Study))+
  labs(x="",y="Coefficient of Variation (%)")+
  theme(legend.position = "right",
        legend.justification = "center",
        legend.key=element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x=element_blank(),
        plot.title = element_text(size=20, hjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin=unit(c(1,-2,-0.3,1), "cm"))

SC_SPIM_ResultsSigma_axLand$Study<-factor(SC_SPIM_ResultsSigma_axLand$Study,
                                      levels=c("Richardson","Algar"))
plot_Sigma_CoV<-ggplot(SC_SPIM_ResultsSigma_axLand[SC_SPIM_ResultsSigma_axLand$DataType!="grpSize1",], 
       aes(x = Model, y = CoV*100,label=Year)) +
  geom_line(aes(group = PairNumber, color=Year),lwd=1.5,alpha=0.5)+
  geom_point(size=4,aes(fill=Year),pch=21, color="black")+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))+
  facet_wrap(vars(Study))+
  labs(x="",y="Coefficient of Variation (%)")+
  theme(strip.text.x = element_blank(),
        legend.position = "right",
        legend.justification = "center",
        legend.key=element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 12),
        plot.title = element_text(size=20, hjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin=unit(c(-0.3,-2,1,1), "cm")) #top right bottom left

#bring in the CoV from the simulations to put alongside
load("simulated_CoV.rda")
sim_long$Study<-"Simulations"
sim_long%>%group_by(Param,Model)%>% summarise(mean=mean(value),
                                                min=min(value),max=max(value))

  
plot_D_CoVb<-ggplot(sim_long[sim_long$Param=="N",],  aes(x = Model, y = value*100,label=Iter)) +
  geom_line(aes(group = Iter),lwd=1.5,alpha=0.2)+#,lty=Year
  geom_point(size=4,pch=20, color="black")+
  # geom_text()+
  scale_y_continuous(limits = c(0, 160),breaks=NULL)+ #breaks = seq(0, 160, by = 20))+#breaks = seq(0, 180, by = 20))+
  facet_wrap(vars(Study))+
  labs(x="",y="Coefficient of Variation (%)")+
  theme(legend.position = "right",
        legend.justification = "center",
        legend.key=element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(size=20, hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        aspect.ratio=18/16,
        plot.margin=unit(c(1,0,-0.3,-0.6), "cm")) #top right bottom left

plot_Sigma_CoVb<-ggplot(sim_long[sim_long$Param=="Sigma",],  aes(x = Model, y = value*100,label=Iter)) +
  geom_line(aes(group = Iter),lwd=1.5,alpha=0.2)+#,lty=Year
  geom_point(size=4,pch=20, color="black")+
  scale_y_continuous(limits = c(0, 100),breaks=NULL)+ # breaks = seq(0, 100, by = 20))+#breaks = seq(0, 180, by = 20))+
  facet_wrap(vars(Study))+
  labs(x="",y="Coefficient of Variation (%)")+
  theme(strip.text.x = element_blank(),
        legend.position = "right",
        legend.justification = "center",
        legend.key=element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y=element_blank(),
        axis.text = element_text(size = 12),
        plot.title = element_text(size=20, hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        aspect.ratio=18/16,
        plot.margin=unit(c(-0.3,0,1,-0.6), "cm")) #top right bottom left

plot_Sigma_CoVb
plot_D_CoVb

plot_AlgarRich_CoV_all<-ggpubr::ggarrange(plot_D_CoV+ rremove("ylab"),
                                          plot_D_CoVb+ rremove("ylab"),
                                          plot_Sigma_CoV+ rremove("ylab"),
                                          plot_Sigma_CoVb+ rremove("ylab"),
                                          common.legend = TRUE,
                                          legend = "bottom")


plot_AlgarRich_CoV_all
annotate_figure(plot_AlgarRich_CoV_all, left = textGrob("Coefficient of variation (CV; %)", rot = 90,
                                                        hjust=0.4, vjust = 2.5, gp = gpar(cex = 1.3)),
                bottom=textGrob("Model",vjust = -5, gp = gpar(cex = 1.3)),
                right = textGrob("Density (D)                     Sigma (\u03c3 ; km)", rot = 270,
                                 hjust=0.7, vjust = 6.5, gp = gpar(cex = 1.3)))





plot_fig2<-plot_AlgarRichv2
plot_fig3<-plot_AlgarRich_CoV_all
plot_figS1<-plot_AlgarRich_allVsgroupSize1v2
