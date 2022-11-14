###############################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------1. CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------2. IMPORT REQUIRED LIBRARIES ------
library(rgdal)
library(raster)
library(coda)
library(nimble)
library(nimbleSCR)
library(stringr)
library(abind)
library(R.utils)
library(sf)
library(fasterize)
library(dplyr)
library(lubridate)
library(stars)
library(RANN)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
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
ResDir <- file.path(thisDir, "results")



all_res <- read.csv(file.path(ResDir, "results_transectsrepetitions.csv"))
def_res <- read.csv(file.path(ResDir, "res5_2.csv"))
def_res["scenario"] <- "full dataset"

all_res <- all_res[,-1]
def_res <- def_res[,-1]

all_res <- rbind(filter(all_res, stat == "means"),filter(all_res, stat == "sd"))
def_res <- rbind(filter(def_res, stat == "means"),filter(def_res, stat == "sd"))


def_colors <- ghibli::ghibli_palette("PonyoMedium")[4]
my_colors <- wes_palette("GrandBudapest2")[1:3]
# my_colors <- c(my_colors, wesanderson::wes_palette("GrandBudapest2")[4:3])

## ------   data wrangling

# means <- filter(all_res, stat == "means")
# sd <- filter(all_res, stat == "sd")
# CV <- filter(all_res, stat == "CV")
# lci <- filter(all_res, stat == "lci")
# uci <- filter(all_res, stat == "uci")


## -----------------------------------------------------------------------------
## ------   N   -----

# a <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5) +
#   geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) + 
#   labs(title="means", fill ="Repetitions simulation", y = "N", x = "Simulation") 
# 


facet_lab <- as_labeller(c('means' = "Means", 
                           'sd' = "Standard Deviation"))


vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  geom_point(data = def_res, size =3, alpha = 1.2, color = def_colors) +
  # scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  labs(title="N - Wolves abundance", fill ="Search Events", y = "N", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal") 
  
ggsave("transectsrepetitionsvl_N_3-6.jpeg", dpi = 300)


bx <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  geom_point(data = def_res, size =3, alpha = 1.2, color = def_colors) +
  # scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  labs(title="N - Wolves abundance", fill ="Search Events", y = "N", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal") 

ggsave("transectsrepetitionsbx_N_3-6.jpeg", dpi = 300)

## -----------------------------------------------------------------------------
## ------   p0   -----

##------   data wrangling


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


x_labs <- c("female RI","female offspring","female other",
  "male RI","male offspring","male other")


## ------   plots

# a <- ggplot(data = p0_means, aes(x=as.character(variable), y=value, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5)  +
#  scale_x_discrete(labels= x_labs) + labs(title="Means", x ="sex-status", y = "p0",fill = "Simulation")
# 

### VIOLIN
vl <- ggplot(data = p0, aes(x=variable, y=value, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2, scale = "width") +
  # geom_boxplot(width = 0.2, color="grey", alpha = 0.2) +
  geom_point(data = p0_52, size=3,alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  labs(title="p0 - baseline detectability", fill ="Search events", y = "p0", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsvl_p0_3-6.jpeg", dpi = 300)

### BOXPLOT

bx <- ggplot(data = p0, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1) +
  geom_point(data = p0_52, size=3, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  labs(title="p0 - baseline detectability", fill ="Search events", y = "p0", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetitionsbx_p0_3-6.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   sigma   -----

## ------   data wrangling



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

vl + geom_violin(alpha = 1.2, scale = "width") +
  geom_point(data = sigma_52, size = 3, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  labs(title="σ - scale parameter", fill ="Search events", y = "sigma", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsvl_sigma_3-6.jpeg", dpi = 300)



bx <- ggplot(data = sigma, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  geom_point(data = sigma_52, size = 3, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y",labeller=facet_lab) + 
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values=my_colors) +
  labs(title="σ - scale parameter", fill ="Search events", y = "sigma", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsbx_sigma_3-6.jpeg", dpi = 300)

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
  geom_point(data = detcov_52, size = 3, alpha = 1.2, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller=facet_lab) + 
  scale_x_discrete(labels= DetCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Search events", y = "Detection Covariates", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsvl_detcov_3-6.jpeg", dpi = 300)



bx <- ggplot(data = detcov, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.2) +
  geom_point(data = detcov_52, size = 3, alpha = 1.2, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller=facet_lab) + 
  scale_x_discrete(labels= DetCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Search events", y = "Detection Covariates", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsbx_detcov_3-6.jpeg", dpi = 300)


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
  geom_point(data = habcov_52, size = 3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller=facet_lab) + 
  scale_x_discrete(labels= HabCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Search events", y = "Habitat Covariates", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsvl_habcov_3-6.jpeg", dpi = 300)



bx <- ggplot(data = habcov, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6,alpha = 1.2) +
  geom_point(data = habcov_52, size = 3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller=facet_lab) + 
  scale_x_discrete(labels= HabCov) +
  scale_fill_manual(values=my_colors) +
  labs(title="Detection Covariates", fill ="Search events", y = "Habitat Covariates", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsbx_habcov_3-6.jpeg", dpi = 300)

## -----------------------------------------------------------------------------
## ------   psi   -----

vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=psi, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2) +
  geom_point(data = def_res,size =3, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab)  +
  scale_fill_manual(values=my_colors) +
  labs(title="ψ - real individual probability", fill ="Search events", y = "psi", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")


ggsave("transectsrepetitionsvl_psi_3-6.jpeg", dpi = 300)


bx <- ggplot(data = all_res, aes(x=as.character(scenario), y=psi, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6,alpha = 1.2) +
  geom_point(data = def_res, size =3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab)  +
  scale_fill_manual(values=my_colors) +
  labs(title="ψ - real individual probability", fill ="Search events", y = "psi", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsbx_psi_3-6.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   rho   -----

vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=rho, fill = as.character(scenario))) 

vl + geom_violin(alpha = 1.2) +
  geom_point(data = def_res, size =3, alpha = 1.2, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) +
  scale_fill_manual(values=my_colors) +
  labs(title="ρ - sex proproportion", fill ="Search events", y = "rho", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsvl_rho_3-6.jpeg", dpi = 300)

bx <- ggplot(data = all_res, aes(x=as.character(scenario), y=rho, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6,alpha = 1.2) +
  geom_point(data = def_res, size =3, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab)  +
  scale_fill_manual(values=my_colors) +
  labs(title="ρ - sex proproportion", fill ="Search events", y = "rho", x = "Number of Search Events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.position="none", legend.box = "horizontal")

ggsave("transectsrepetitionsbx_rho_3-6.jpeg", dpi = 300)

