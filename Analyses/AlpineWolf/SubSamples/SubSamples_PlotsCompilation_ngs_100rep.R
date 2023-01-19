###############################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------1. CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------2. IMPORT REQUIRED LIBRARIES ------

library(stringr)
library(abind)
library(R.utils)
library(dplyr)
library(tidyr)
library(purrr)
library(reshape2)
library(lubridate)
library(stars)
library(RANN)
library(gridExtra)
library(MetBrewer)
library(fs)
library(ggplot2)
library(ggpmisc)
library(wesanderson)
library(ghibli)
library(hrbrthemes)


## ------ 3. SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ 4.  SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))


## -----------------------------------------------------------------------------
## ------ 5. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.SubSample.ResultsOutput"
thisDir <- file.path(analysisDir, modelName)


## -----------------------------------------------------------------------------
## ------   6. PLOT RESULTS -----

ResDir <- file.path(thisDir, "results")


# Load results
all_res <- read.csv(file.path(ResDir, "results_transectsubsample_25-50-75.csv"))
def_res <- read.csv(file.path(ResDir, "res5_2.csv"))
def_res["scenario"] <- "full dataset"

# Leave out first column

all_res <- all_res[,-1]
def_res <- def_res[,-1]

# Just keep means and sd

# all_res <- rbind(filter(all_res, stat == "means"),filter(all_res, stat == "sd"))
# def_res <- rbind(filter(def_res, stat == "means"),filter(def_res, stat == "sd"))

# Set colours

def_colors <- ghibli::ghibli_palette("PonyoMedium")[4]
my_colors <- wes_palette("GrandBudapest2")[1:2]
my_colors <- c(my_colors, wesanderson::wes_palette("GrandBudapest2")[4:3])



## -----------------------------------------------------------------------------
## ------   N   -----

# Reordering stat factor levels
def_res$stat <- factor(def_res$stat, levels = c("means", "sd","CV","lci", "uci"))
all_res$stat <- factor(all_res$stat, levels = c("means", "sd","CV","lci", "uci"))


# Create nice labels for the stats in the plots
facet_lab <- as_labeller(c('means' = "Means", 
                           'sd' = "Standard Deviation",
                           'CV' = "Coefficient of Variation",
                           'lci'="Lower Credibility Interval",
                           'uci'="Upper Credibility Interval"))

x_labs <- c("25%", "50%", "75%", "100%")


# VL violin plots

# vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2) +
#   facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
#   geom_point(data = def_res, size =4, alpha = 1.2, color = def_colors) +
#   # scale_x_discrete(labels= x_labs) +
#   scale_fill_manual(values=my_colors) +
#   labs(title="N - Wolves abundance", fill ="NGS samples %", y = "N", x = "NGS samples %") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size=11),
#         axis.text.x = element_text(size=11,hjust=1),
#         axis.text.y = element_text(size=11),
#         legend.position="none", legend.box = "horizontal") 

# uncomment to save
# ggsave("ngsvl_N_25-50-75.jpeg", dpi = 300)

# BX box plots

bx <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  geom_jitter(alpha = 0.2, color = "grey") + 
  geom_point(data = def_res, size =4, alpha = 1.2, color = def_colors) +
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  # expand_limits( y = 0) +

  labs(title="N - Wolves abundance", fill ="NGS samples %", y = "N", x = "NGS samples") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="bottom", legend.box = "horizontal") 

# ggsave("ngsbx_N_25-50-75.jpeg", dpi = 300)

## -----------------------------------------------------------------------------
## ------   p0   -----

##------   Data wrangling


p0 <- all_res[, grep("p0", colnames(all_res))]
p0["stat"] <- all_res$stat
p0["scenario"] <- all_res$scenario
p0 <- melt(p0,id.vars=c("scenario","stat"))
# p01 <- p0 %>% 
#   filter(stat != "lci"| stat != "uci")


p0_52 <- def_res[, grep("p0", colnames(def_res))]
p0_52["stat"] <- def_res$stat
p0_52["scenario"] <- def_res$scenario
p0_52 <- melt(p0_52,id.vars=c("scenario","stat"))

# Labs for sexes and statuses
x_labs_t <- c("F RI","F offspring","F other",
              "M RI","M offspring","M other")


## ------   Plots


### VIOLIN
# vl <- ggplot(data = p0, aes(x=variable, y=value, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2) +
#   # geom_boxplot(width = 0.2, color="grey", alpha = 0.2) +
#   geom_point(data = p0_52, size =3,alpha = 1.2,color = def_colors) +
#   facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
#   scale_x_discrete(labels= x_labs_t) +
#   scale_fill_manual(values=my_colors) +
#   labs(title="p0 - baseline detectability", fill ="NGS samples %", y = "p0", x = "NGS samples %") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size=11),
#         axis.text.x = element_text(size=11,angle=45,hjust=1),
#         axis.text.y = element_text(size=11),
#         legend.position="none", legend.box = "horizontal")

# ggsave("ngsvl_p0_25-50-75.jpeg", dpi = 300)


### BOXPLOT

bx <- ggplot(data = p0, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1) +
  geom_boxplot(data = p0_52, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  scale_x_discrete(labels= x_labs_t) +
  scale_fill_manual(values=my_colors) +
  labs(title="p0 - baseline detectability", fill ="NGS samples %", y = "p0", x = "NGS samples %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=10),
        axis.text.x = element_text(size=11,angle=45,hjust=1),
        axis.text.y = element_text(size=10),
        legend.position="right", legend.box = "vertical")

# ggsave("ngsbx_p0_25-50-75.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   sigma   -----

## ------   Data wrangling



sigma <- all_res[, grep("sigma", colnames(all_res))]
sigma["stat"] <- all_res$stat
sigma["scenario"] <- all_res$scenario
sigma <- melt(sigma,id.vars=c("scenario","stat"))

sigma_52 <- def_res[, grep("sigma", colnames(def_res))]
sigma_52["stat"] <- def_res$stat
sigma_52["scenario"] <- def_res$scenario
sigma_52 <- melt(sigma_52,id.vars=c("scenario","stat"))

## ------   plots

# a <- ggplot(data = sigma_means, aes(x=as.character(variable), y=value, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5) +
#   scale_x_discrete(labels= x_labs) +labs(title="Means", x ="sex-status", y = "sigma", fill = "Simulation")


# vl <- ggplot(data = sigma, aes(x=variable, y=value, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2) +
#   #  add scale = "width" in geom_violin
#   geom_boxplot(data = sigma_52, alpha = 1.2,color = def_colors) +
#   facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
#   scale_x_discrete(labels= x_labs) +
#   scale_fill_manual(values=my_colors) +
#   labs(title="σ - scale parameter", fill ="NGS samples %", y = "sigma", x = "NGS samples %") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size=11),
#         axis.text.x = element_text(size=11,angle=45,hjust=1),
#         axis.text.y = element_text(size=11),
#         legend.position="none", legend.box = "horizontal")

# ggsave("ngsvl_sigma_25-50-75.jpeg", dpi = 300)



bx <- ggplot(data = sigma, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  geom_boxplot(data = sigma_52, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
  scale_x_discrete(labels= x_labs_t) +
  scale_fill_manual(values=my_colors) +
  labs(title="σ - scale parameter", fill ="NGS samples %", y = "sigma", x = "NGS samples %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,angle=45,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="right", legend.box = "vertical")

# ggsave("ngsbx_sigma_25-50-75.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   DetCov   -----

## ------   data wrangling

DetCov <- c("Transects Length", "Operator Experience", "Snow", "East/West", "Log Human Pop")

detcov <- all_res[,grep("betaDet", colnames(all_res))]
detcov["stat"] <- all_res$stat
detcov["scenario"] <- all_res$scenario
detcov <- melt(detcov,id.vars=c("scenario","stat"))

detcov_52 <- def_res[,grep("betaDet", colnames(def_res))]
detcov_52["stat"] <- def_res$stat
detcov_52["scenario"] <- def_res$scenario
detcov_52 <- melt(detcov_52,id.vars=c("scenario","stat"))

## ------   plots

# a <- ggplot(data = detcov_means, aes(x=as.character(variable), y=value, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5) +
#   scale_x_discrete(labels= DetCov) +labs(title="Means", x ="sex-status", y = "Detection Covariates", fill = "Simulation")

# 
# vl <- ggplot(data = detcov, aes(x=variable, y=value, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2) +
#   # geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
#   geom_point(data = detcov_52, size =3, alpha = 1.2, color = def_colors) +
#   facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
#   scale_x_discrete(labels= DetCov) +
#   scale_fill_manual(values=my_colors) +
#   labs(title="Detection Covariates", fill ="Search events", y = "Detection Covariates", x = "NGS samples %") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size=11),
#         axis.text.x = element_text(size=11,angle=45,hjust=1),
#         axis.text.y = element_text(size=11),
#         legend.position="none", legend.box = "horizontal")

# ggsave("ngsvl_detcov_25-50-75.jpeg", dpi = 300)

bx <- ggplot(data = detcov, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.2) +
  geom_point(data = detcov_52, size =3, alpha = 1.2, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
  scale_x_discrete(labels= DetCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Search events", y = "Detection Covariates", x = "NGS samples %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

# ggsave("ngsbx_detcov_25-50-75.jpeg", dpi = 300)

## -----------------------------------------------------------------------------
## ------   HabCov   -----

## ------   data wrangling
HabCov <- c("Bare Rocks","Herbaceous", "Forest","Human Pop","Wolf presence")

habcov <- all_res[,grep("betaDet", colnames(all_res))]
habcov["stat"] <- all_res$stat
habcov["scenario"] <- all_res$scenario
habcov <- melt(habcov,id.vars=c("scenario","stat"))

habcov_52 <- def_res[,grep("betaDet", colnames(def_res))]
habcov_52["stat"] <- def_res$stat
habcov_52["scenario"] <- def_res$scenario
habcov_52 <- melt(habcov_52,id.vars=c("scenario","stat"))




# vl <- ggplot(data = habcov, aes(x=variable, y=value, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2) +
#   # geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
#   geom_point(data = habcov_52, size =3, alpha = 1.2, color=def_colors) +
#   facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
#   scale_x_discrete(labels= HabCov) +
#   scale_fill_manual(values=my_colors) +
#   labs(title="Detection Covariates", fill ="Search events", y = "Habitat Covariates", x = "NGS samples %") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size=11),
#         axis.text.x = element_text(size=11,angle=45,hjust=1),
#         axis.text.y = element_text(size=11),
#         legend.position="none", legend.box = "horizontal")

# ggsave("ngsvl_habcov_25-50-75.jpeg", dpi = 300)



bx <- ggplot(data = habcov, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6,alpha = 1.2) +
  geom_point(data = habcov_52, size =3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  scale_x_discrete(labels= HabCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Search events", y = "Habitat Covariates", x = "NGS samples %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,angle=45,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")

# ggsave("ngsbx_habcov_25-50-75.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   psi   -----

# vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=psi, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2) +
#   geom_point(data = def_res,size =3, alpha = 1.2,color = def_colors) +
#   facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab)  +
#   scale_fill_manual(values=my_colors) +
#   labs(title="ψ - real individual probability", fill ="Search events", y = "psi", x = "NGS samples %") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size=11),
#         axis.text.x = element_text(size=11,hjust=1),
#         axis.text.y = element_text(size=11),
#         legend.position="none", legend.box = "horizontal")


# ggsave("ngsvl_psi_25-50-75.jpeg", dpi = 300)


bx <- ggplot(data = all_res, aes(x=as.character(scenario), y=psi, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6,alpha = 1.2) +
  geom_point(data = def_res, size =3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab)  +
  scale_fill_manual(values=my_colors) +
  labs(title="ψ - real individual probability", fill ="Search events", y = "psi", x = "NGS samples %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")

# ggsave("ngsbx_habcov_25-50-75.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   rho   -----

# vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=rho, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2) +
#   geom_point(data = def_res, size =3, alpha = 1.2,color=def_colors) +
#   facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) +
#   scale_fill_manual(values=my_colors) +
#   labs(title="ρ - sex proproportion", fill ="Search events", y = "rho", x = "NGS samples %") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size=11),
#         axis.text.x = element_text(size=11,hjust=1),
#         axis.text.y = element_text(size=11),
#         legend.position="none", legend.box = "horizontal")

# ggsave("ngsvl_rho_25-50-75.jpeg", dpi = 300)

bx <- ggplot(data = all_res, aes(x=as.character(scenario), y=rho, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6,alpha = 1.2) +
  geom_point(data = def_res, size =3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab)  +
  scale_fill_manual(values=my_colors) +
  labs(title="ρ - sex proproportion", fill ="Search events", y = "rho", x = "NGS samples %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")



# ggsave("ngsbx_rho_25-50-75.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   7. PLOT RESULTS II -----


facet_lab_ii <- as_labeller(c("25"= "25%", "50"="50%", "75"="75%"))


# IDs
all_in <- read.csv(file.path(ResDir, "transects_input.csv"))
all_in <- subset(all_in, variable == "IDs")

all_in <- all_in[,-1]

colnames(all_in) <- c('N','scenario')
all_in["source"] <- "inputs"


# all_in_25 <- subset(all_in_1, scenario == "25")
# all_in_50 <- subset(all_in_1, scenario == "50")
# all_in_75 <- subset(all_in_1, scenario == "75")


all_res <- subset(all_res, stat == "means")
all_res <- select(all_res, N, scenario)
all_res["source"] <- "outputs"

# all_res_N_25 <- subset(all_res_N_tot, V2 == "25")
# all_res_N_50 <- subset(all_res_N_tot, V2 == "50")
# all_res_N_75 <- subset(all_res_N_tot, V2 == "75")

dataplot <- rbind(all_res,all_in)
dataplot <- mutate(dataplot, ID = row_number())
dataplot.1 <-  dcast(dataplot, ID + scenario ~ source, value.var =  "N", fun.aggregate = sum)
# dataplot.2 <- reshape(dataplot, idvar = "scenario", timevar = "source", direction = "wide")

# write.csv(dataplot.1, "in_res_TRNS.csv")
dataplot <- read.csv(file.path(ResDir, "in_res_NGS.csv"))


bx <- ggplot(data = dataplot, aes(x=outputs, y=inputs)) 
bx +   facet_grid(~scenario)  +
  geom_jitter(alpha = 0.6)  +
  stat_poly_line() +
  stat_poly_eq() +
  labs(title="NGS Subsampling", y = "Wolf genotypes", x = "Wolves abundace estimates") +
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11))
ggsave("IDs_nest_NGS.jpeg", dpi = 300)

## -----------------------------------------------------------------------------
# Detections
  
all_in <- read.csv(file.path(ResDir, "transects_input.csv"))
all_in <- subset(all_in, variable == "NDetections")

all_in <- all_in[,-1]

colnames(all_in) <- c('N','scenario')
all_in["source"] <- "inputs"



all_res <- subset(all_res, stat == "means")
all_res <- select(all_res, N, scenario)
all_res["source"] <- "outputs"



dataplot <- rbind(all_res,all_in)
dataplot <- mutate(dataplot, ID = row_number())
dataplot.1 <-  dcast(dataplot, ID + scenario ~ source, value.var =  "N", fun.aggregate = sum)
# dataplot.2 <- reshape(dataplot, idvar = "scenario", timevar = "source", direction = "wide")

# write.csv(dataplot.1, "in_det_res_TRNS.csv")

dataplot <- read.csv(file.path(ResDir, "in_det_res_NGS.csv"))


bx <- ggplot(data = dataplot, aes(x=outputs, y=inputs)) 
bx +   facet_grid(~scenario)  +
  geom_jitter(alpha = 0.6)  +
  stat_poly_line() +
  stat_poly_eq() +
  labs(title="NGS Subsampling", y = "Wolf detections", x = "Wolf abundace estimates") +
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11))

ggsave("det_nest_NGS.jpeg", dpi = 300)

