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
library(lubridate)
library(stars)
library(RANN)
library(gridExtra)
library(MetBrewer)
library(fs)
library(purrr)
library(ggplot2)
library(reshape2)
library(ghibli)
library(wesanderson)
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
# Load results csv
ResDir <- file.path(thisDir, "results")

all_res <- read.csv(file.path(thisDir, "results/results_transectsubsample_prov_tr.csv"))
def_res <- read.csv(file.path(thisDir, "results/res5_2.csv"))
def_res["scenario"] <- "total dataset"

all_res <- all_res[,-1]
def_res <- def_res[,-1]

# Just keep means and sd

all_res_m <- rbind(filter(all_res, stat == "means"))
all_res_sd <- rbind(filter(all_res, stat == "sd"))

def_res_m <- rbind(filter(def_res, stat == "means"))
def_res_sd <- rbind(filter(def_res, stat == "sd"))

# Set colours
def_colors <- ghibli::ghibli_palette("PonyoMedium")[4]
my_colors <- wes_palette("GrandBudapest2")[1:2]
my_colors <- c(my_colors, wesanderson::wes_palette("GrandBudapest2")[4:3])


## -----------------------------------------------------------------------------
## ------   N   -----

# Create nice labels for the stats in the plots
facet_lab <- as_labeller(c('means' = "Means", 
                           'sd' = "Standard Deviation"))


x_labs <- c("25%", "50%", "75%", "100%")
# VL violin plots

# vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2) +
#   facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
#   geom_point(data = def_res, size =4, alpha = 1.2, color = def_colors) +
#   scale_x_discrete(labels= x_labs) +
#   scale_fill_manual(values=my_colors) +
#   labs(title="N - Wolves abundance", fill ="Transects %", y = "N", x = "Transects %") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size=11),
#         axis.text.x = element_text(size=11,hjust=1),
#         axis.text.y = element_text(size=11),
#         legend.position="none", legend.box = "horizontal") 

# uncomment to save
# ggsave("transvl_N_25-50-75.jpeg", dpi = 300)

# BX box plots

bx <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  geom_point(data = def_res, size =4, alpha = 1.2, color = def_colors) +
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  labs(title="N - Wolves abundance", fill ="Transectss %", y = "N", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal") 


# ----- 3 % 6 repetitions transects MEANS

bx <- ggplot(data = all_res_m, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  facet_grid(~ rep)  +
  geom_point(data = def_res_m, size = 4, alpha = 1.2, color = def_colors) +
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  labs(title="N - Wolves abundance", fill ="Transectss %", y = "N", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal") 
# ggsave("transbx_N_25-50-75.jpeg", dpi = 300)

# ----- 3 % 6 repetitions transects SD


bx <- ggplot(data = all_res_m, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  facet_grid(~ rep)  +
  geom_point(data = def_res_m, size = 4, alpha = 1.2, color = def_colors) +
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  labs(title="N - Wolves abundance", fill ="Transectss %", y = "N", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal") 

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
vl <- ggplot(data = p0, aes(x=variable, y=value, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2) +
  # geom_boxplot(width = 0.2, color="grey", alpha = 0.2) +
  geom_point(data = p0_52, size =3,alpha = 1,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  scale_x_discrete(labels= x_labs_t) +
  scale_fill_manual(values=my_colors) +
  labs(title="p0 - baseline detectability", fill ="Transects %", y = "p0", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,angle=45,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")

# ggsave("transvl_p0_25-50-75.jpeg", dpi = 300)


### BOXPLOT

bx <- ggplot(data = p0, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  geom_boxplot(data = p0_52, width = 0.6, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  scale_x_discrete(labels= x_labs_t) +
  scale_fill_manual(values=my_colors) +
  labs(title="p0 - baseline detectability", fill ="Transects %", y = "p0", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=10,angle=45,hjust=1),
        axis.text.y = element_text(size=10),
        legend.position="none", legend.box = "vertical")

ggsave("transbx_prov_p0_25-50-75.jpeg", dpi = 300)


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


vl <- ggplot(data = sigma, aes(x=variable, y=value, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2,scale = "width") +
  #  add scale = "width" in geom_violin +
  facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
  scale_x_discrete(labels= x_labs_t) +
  scale_fill_manual(values=my_colors) +
  geom_boxplot(data = sigma_52, color = def_colors) +
  labs(title="σ - scale parameter", fill ="Transects %", y = "sigma", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,angle=45,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")

# ggsave("transvl_sigma_25-50-75.jpeg", dpi = 300)



bx <- ggplot(data = sigma, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  geom_boxplot(data = sigma_52, width = 0.6, alpha = 0.5, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
  scale_x_discrete(labels= x_labs_t) +
  scale_fill_manual(values=my_colors) +
  labs(title="σ - scale parameter", fill ="Transects %", y = "sigma", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=10),
        axis.text.x = element_text(size=11,angle=45,hjust=1),
        axis.text.y = element_text(size=10),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transbx_prov_sigma_25-50-75.jpeg", dpi = 300)


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


vl <- ggplot(data = detcov, aes(x=variable, y=value, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2) +
  # geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = detcov_52, size =3, alpha = 1.2, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
  scale_x_discrete(labels= DetCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Transects %", y = "Detection Covariates", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,angle=45,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")

ggsave("transvl_detcov_25-50-75.jpeg", dpi = 300)



bx <- ggplot(data = detcov, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.2) +
  geom_point(data = detcov_52, size =3, alpha = 1.2, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
  scale_x_discrete(labels= DetCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Transects %", y = "Detection Covariates", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="none", legend.box = "horizontal")

ggsave("transbx_detcov_25-50-75.jpeg", dpi = 300)

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




vl <- ggplot(data = habcov, aes(x=variable, y=value, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2) +
  # geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = habcov_52, size =3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  scale_x_discrete(labels= HabCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Transects %", y = "Habitat Covariates", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,angle=45,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")

ggsave("transvl_habcov_25-50-75.jpeg", dpi = 300)



bx <- ggplot(data = habcov, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6,alpha = 1.2) +
  geom_point(data = habcov_52, size =3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  scale_x_discrete(labels= HabCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Transects %", y = "Habitat Covariates", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,angle=45,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")

ggsave("transbx_habcov_25-50-75.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   psi   -----

vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=psi, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2) +
  geom_point(data = def_res,size =3, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab)  +
  scale_fill_manual(values=my_colors) +
  labs(title="ψ - real individual probability", fill ="Transects %", y = "psi", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")


ggsave("transvl_psi_25-50-75.jpeg", dpi = 300)


bx <- ggplot(data = all_res, aes(x=as.character(scenario), y=psi, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6,alpha = 1.2) +
  geom_point(data = def_res, size =3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab)  +
  scale_fill_manual(values=my_colors) +
  labs(title="ψ - real individual probability", fill ="Transects %", y = "psi", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")

ggsave("transbx_psi_25-50-75.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   rho   -----

vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=rho, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2) +
  geom_point(data = def_res, size =3, alpha = 1.2,color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) +
  scale_fill_manual(values=my_colors) +
  labs(title="ρ - sex proproportion", fill ="Transects %", y = "rho", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")

ggsave("transvl_rho_25-50-75.jpeg", dpi = 300)

bx <- ggplot(data = all_res, aes(x=as.character(scenario), y=rho, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6,alpha = 1.2) +
  geom_point(data = def_res, size =3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab)  +
  scale_fill_manual(values=my_colors) +
  labs(title="ρ - sex proproportion", fill ="Transects %", y = "rho", x = "Transects %") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=11,hjust=1),
        axis.text.y = element_text(size=11),
        legend.position="none", legend.box = "horizontal")



ggsave("transbx_rho_25-50-75.jpeg", dpi = 300)

