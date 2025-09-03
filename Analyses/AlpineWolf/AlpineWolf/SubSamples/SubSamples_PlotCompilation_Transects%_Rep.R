###############################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------1. CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------2. IMPORT REQUIRED LIBRARIES ------
# library(rgdal)
# library(raster)
# library(coda)
# library(nimble)
# library(nimbleSCR)
library(stringr)
library(abind)
library(R.utils)
library(sf)
library(fasterize)
library(dplyr)
# library(lubridate)
# library(stars)
# library(RANN)
# library(Rcpp)
# library(RcppArmadillo)
# library(RcppProgress)
library(gridExtra)
library(MetBrewer)
library(fs)
library(purrr)
library(ggplot2)
library(reshape2)
library(ghibli)
library(wesanderson)
library(hrbrthemes)
library(RColorBrewer)
library(cowplot)
library(forcats)


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

all_res_t <- read.csv(file.path(ResDir, "results_transectsubsample_prov.csv"))
all_res_t["subsampling"] <- "Transects Spatial"
all_res_t["rep"] <- "All Transects repetitions"

all_res_tr <- read.csv(file.path(ResDir, "results_transectsrepetitions.csv"))
all_res_tr["subsampling"] <- "Transects Spatial and Intensity"
colnames(all_res_tr)[34] <- "rep"
all_res_tr["scenario"] <- "100"
all_res_tr$rep[all_res_tr$rep=="3"]<-"3 repetitions"
all_res_tr$rep[all_res_tr$rep=="6"]<-"6 repetitions"


all_res_t_tr <- read.csv(file.path(ResDir, "results_transectsubsample_prov_tr.csv"))
all_res_t_tr["subsampling"] <- "Transects Spatial and Intensity"
all_res_t_tr$rep[all_res_t_tr$rep=="3"]<-"3 repetitions"
all_res_t_tr$rep[all_res_t_tr$rep=="6"]<-"6 repetitions"


def_res <- read.csv(file.path(ResDir, "res5_2.csv"))
def_res["scenario"] <- "full dataset"
def_res["rep"] <- "full dataset"

all_res <- rbind(all_res_t, all_res_tr,all_res_t_tr)


all_res <- all_res[,-1]
def_res <- def_res[,-1]

# all_res <- rbind(filter(all_res, stat == "means"),filter(all_res, stat == "sd"))
# def_res <- rbind(filter(def_res, stat == "means"),filter(def_res, stat == "sd"))

# Set colors
def_colors <- ghibli::ghibli_palette("PonyoMedium")[4]
my_colors <- wes_palette("GrandBudapest2")[1:2]
my_colors1 <- wes_palette("GrandBudapest2")[4]

my_colors2 <- ghibli::ghibli_palette("PonyoMedium")[3:4]

my_colors <- c(my_colors, my_colors1, my_colors2)

## ------   data wrangling

means <- filter(all_res, stat == "means")
sd <- filter(all_res, stat == "sd")
CV <- filter(all_res, stat == "CV")
lci <- filter(all_res, stat == "lci")
uci <- filter(all_res, stat == "uci")

def_res_m <- filter(def_res, stat == "means")
def_res_sd <- filter(def_res, stat == "sd")
def_res_cv <- filter(def_res, stat == "CV")
def_res_lci <- filter(def_res, stat == "lci")
def_res_uci <- filter(def_res, stat == "uci")


## -----------------------------------------------------------------------------
## ------   N   -----

# a <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5) +
#   geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) + 
#   labs(title="means", fill ="Repetitions simulation", y = "N", x = "Simulation") 
# 
# means$scenario <- replace(means$scenario, means$scenario == 25, '25 %')
# means$scenario <- replace(means$scenario, means$scenario == 50, '50 %')
# means$scenario <- replace(means$scenario, means$scenario == 75, '75 %')
# means$scenario <- replace(means$scenario, means$scenario == 100, '100 %')


# means$rep <- replace(means$rep, means$rep == 3, '3 search events')
# means$rep <- replace(means$rep, means$rep == 6, '6 search events')

# 
#  vl <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario)))
# 
#  vl + geom_violin(alpha = 1.2) +
#    facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) +
#    geom_point(data = def_res, size =3, alpha = 1.2, color = def_colors) +
#    # scale_x_discrete(labels= x_labs) +
#    scale_fill_manual(values=my_colors) +
#    labs(title="N - Wolves abundance", fill ="Search Events", y = "N", x = "Number of Search Events") +
#    theme_ipsum() +
#    theme(plot.title = element_text(size = 11),
#         axis.text.x = element_text(size = 11, hjust = 1),
#         axis.text.y = element_text(size = 11),
#         legend.position="none", legend.box = "horizontal")
# 
#  vl 
# ggsave("transectsrepetitionsvl_N_3-6.jpeg", dpi = 300)


# data_hline <- data.frame(group = unique(def_res_m$scenario),  # Create data for lines
#                          hline = def_res_m)
# data_hline_lci <- data.frame(group = unique(def_res_lci$scenario),  # Create data for lines
#                          hline = def_res_lci)
# data_hline_uci <- data.frame(group = unique(def_res_lci$scenario),  # Create data for lines
#                          hline = def_res_uci)
# 

# ------ AV$SD -----

bx <- means %>%
  mutate(scenario = fct_relevel(scenario,"25", "50", "75", "100")) %>%
  ggplot(aes(x=scenario, y=N, fill = scenario)) +
  geom_violin(width = 0.6, alpha = 1.5) +
  geom_hline(yintercept=def_res_m$N, linetype="dashed", color = "red") +
  facet_grid(vars(rep)) + 
  # geom_hline(data = data_hline,
  #            aes(yintercept = hline.N),
  #                size = 0.2) +
  # geom_hline(data = data_hline_lci,
  #            aes(yintercept = hline.N), 
  #                linetype = "dashed", 
  #                size=0.2) + 
  # geom_hline(data = data_hline_uci,
  #            aes(yintercept = hline.N),
  #            linetype = "dashed", 
  #            size=0.2) +
  # geom_point(data = def_res_m, size =4, alpha = 1.2, color = def_colors) +
  # scale_x_discrete(labels= x_labs) +
  stat_summary(fun = median, geom = "point", color = "white",
             position = position_dodge(0.9)) +
  scale_fill_manual(values=my_colors, name = "Transects Spatial Subsampling",labels=c('25%', '50%', '75%', "100%")) +
  labs(y = "N", x = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=0),
        axis.text.y = element_text(size=10))



bx_sd <- sd %>%
  mutate(scenario = fct_relevel(scenario,"25", "50", "75", "100")) %>%
  ggplot(aes(x=scenario, y=N, fill = scenario)) +
  geom_violin(width = 0.6, alpha = 1.5) +
  geom_hline(yintercept=def_res_sd$N, linetype="dashed", color = "red") +
  facet_grid(vars(rep)) + 
  # geom_hline(data = data_hline,
  #            aes(yintercept = hline.N),
  #                size = 0.2) +
  # geom_hline(data = data_hline_lci,
  #            aes(yintercept = hline.N), 
  #                linetype = "dashed", 
  #                size=0.2) + 
  # geom_hline(data = data_hline_uci,
  #            aes(yintercept = hline.N),
  #            linetype = "dashed", 
  #            size=0.2) +
# geom_point(data = def_res_m, size =4, alpha = 1.2, color = def_colors) +
# scale_x_discrete(labels= x_labs) +
stat_summary(fun = median, geom = "point", color = "white",
             position = position_dodge(0.9)) +
  scale_fill_manual(values=my_colors, name = "Transects Spatial Subsampling",labels=c('25%', '50%', '75%', "100%")) +
  labs(y = "Standard Deviation", x = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=0),
        axis.text.y = element_text(size=9))


plot_grid(bx, bx_sd, labels = "AUTO")

# ------ AV$CV -----

bx <- means %>%
  mutate(scenario = fct_relevel(scenario,"25", "50", "75", "100")) %>%
  ggplot(aes(x=scenario, y=N, fill = scenario)) +
  geom_violin(width = 0.6, alpha = 1.5) +
  geom_hline(yintercept=def_res_m$N, linetype="dashed", color = "red") +
  facet_grid(vars(rep)) + 
  # geom_hline(data = data_hline,
  #            aes(yintercept = hline.N),
  #                size = 0.2) +
  # geom_hline(data = data_hline_lci,
  #            aes(yintercept = hline.N), 
  #                linetype = "dashed", 
  #                size=0.2) + 
  # geom_hline(data = data_hline_uci,
  #            aes(yintercept = hline.N),
  #            linetype = "dashed", 
  #            size=0.2) +
# geom_point(data = def_res_m, size =4, alpha = 1.2, color = def_colors) +
# scale_x_discrete(labels= x_labs) +
stat_summary(fun = median, geom = "point", color = "white",
             position = position_dodge(0.9)) +
  scale_fill_manual(values=my_colors, name = "Transects Spatial Subsampling",labels=c('25%', '50%', '75%', "100%")) +
  labs(y = "N", x = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=0),
        axis.text.y = element_text(size=10))



bx_CV <- CV %>%
  mutate(scenario = fct_relevel(scenario,"25", "50", "75", "100")) %>%
  ggplot(aes(x=scenario, y=N, fill = scenario)) +
  geom_violin(width = 0.6, alpha = 1.5) +
  geom_hline(yintercept=def_res_cv$N, linetype="dashed", color = "red") +
  facet_grid(vars(rep)) + 
  # geom_hline(data = data_hline,
  #            aes(yintercept = hline.N),
  #                size = 0.2) +
  # geom_hline(data = data_hline_lci,
  #            aes(yintercept = hline.N), 
  #                linetype = "dashed", 
  #                size=0.2) + 
  # geom_hline(data = data_hline_uci,
  #            aes(yintercept = hline.N),
  #            linetype = "dashed", 
  #            size=0.2) +
# geom_point(data = def_res_m, size =4, alpha = 1.2, color = def_colors) +
# scale_x_discrete(labels= x_labs) +
stat_summary(fun = median, geom = "point", color = "white",
             position = position_dodge(0.9)) +
  scale_fill_manual(values=my_colors, name = "Transects Spatial Subsampling",labels=c('25%', '50%', '75%', "100%")) +
  labs(y = "CV", x = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=0),
        axis.text.y = element_text(size=9))


plot_grid(bx, bx_CV, labels = "AUTO")

# ggsave("transects-and-repetitionsbx.jpeg", dpi = 300)

## -----------------------------------------------------------------------------
## ------   p0   -----

##------   data wrangling


p0 <- all_res[, grep("p0", colnames(all_res))]
p0["stat"] <- all_res$stat
p0["scenario"] <- all_res$scenario
p0["rep"] <- all_res$rep
p0["subsampling"] <- all_res$subsampling

p0 <- reshape2::melt(p0,id.vars=c("scenario","stat","rep","subsampling"))
# p01 <- p0 %>% 
#   filter(stat != "lci"| stat != "uci")


p0_52 <- def_res[, grep("p0", colnames(def_res))]
p0_52["stat"] <- def_res$stat
p0_52["scenario"] <- def_res$scenario
p0_52["rep"] <- def_res$rep
p0_52["subsampling"] <- def_res$subsampling
pp0_52 <- reshape2::melt(p0_52,id.vars=c("scenario","stat"))


x_labs_t <- c("F RI","F offspring","F other",
              "M RI","M offspring","M other")


## ------   filter

p0means <- filter(p0, stat == "means")
p0sd <- filter(p0, stat == "sd")
p0CV <- filter(p0, stat == "CV")
p0lci <- filter(p0, stat == "lci")
p0uci <- filter(p0, stat == "uci")

p0_52_m <- filter(pp0_52, stat == "means")
p0_52_sd <- filter(pp0_52, stat == "sd")
p0_52_cv <- filter(pp0_52, stat == "CV")
p0_52_lci <- filter(pp0_52, stat == "lci")
p0_52_uci <- filter(pp0_52, stat == "uci")

## ------   plots

# a <- ggplot(data = p0_means, aes(x=as.character(variable), y=value, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5)  +
#  scale_x_discrete(labels= x_labs) + labs(title="Means", x ="sex-status", y = "p0",fill = "Simulation")
# 

### VIOLIN
# vl <- ggplot(data = p0, aes(x=variable, y=value, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2, scale = "width") +
#   # geom_boxplot(width = 0.2, color="grey", alpha = 0.2) +
#   geom_point(data = p0_52, size=3,alpha = 1.2,color = def_colors) +
#   facet_wrap(vars(stat), scales = "free_y", labeller = facet_lab) + 
#   scale_x_discrete(labels= x_labs_t) +
#   scale_fill_manual(values=my_colors) +
#   labs(title="p0 - baseline detectability", fill ="Search events", y = "p0", x = "Number of Search Events") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size = 11),
#         axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
#         axis.text.y = element_text(size = 11),
#         legend.position="none", legend.box = "horizontal")
# 
# ggsave("transectsrepetitionsvl_p0_3-6.jpeg", dpi = 300)

### BOXPLOT

# ------ AV$CV -----

p0m <- p0means %>%
  mutate(scenario = fct_relevel(scenario,"25", "50", "75", "100")) %>%
  ggplot(aes(x=variable, y=value, fill = scenario)) +
  geom_violin(width = 0.6, alpha = 1.5) +
  geom_hline(yintercept=pp0_52$value, linetype="dashed", color = "red") +
  facet_grid(vars(rep)) + 
  stat_summary(fun = median, geom = "point", color = "white",
             position = position_dodge(0.9)) +
  scale_fill_manual(values=my_colors, name = "Transects Spatial Subsampling",labels=c('25%', '50%', '75%', "100%")) +
  labs(y = "N", x = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=0),
        axis.text.y = element_text(size=10))


p0_CV <- p0CV %>%
  mutate(scenario = fct_relevel(scenario,"25", "50", "75", "100")) %>%
  ggplot(aes(x=variable, y=value, fill = scenario)) +
  geom_violin(width = 0.6, alpha = 1.5) +
  geom_hline(yintercept=pp0_52_cv$value, linetype="dashed", color = "red") +
  facet_grid(vars(rep)) + 
  stat_summary(fun = median, geom = "point", color = "white",
             position = position_dodge(0.9)) +
  scale_fill_manual(values=my_colors, name = "Transects Spatial Subsampling",labels=c('25%', '50%', '75%', "100%")) +
  labs(y = "CV", x = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=0),
        axis.text.y = element_text(size=9))


plot_grid(p0m, p0_CV, labels = "AUTO")







bx <- ggplot(data = p0$stat, aes(x=variable, y=value, fill = as.character(scenario)))

bx + geom_boxplot(width = 0.6, alpha = 1) +
  geom_boxplot(data = p0_52, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y") + 
  scale_x_discrete(labels= x_labs_t) +
  scale_fill_manual(values=my_colors) +
  labs(title="p0 - baseline detectability", fill ="Maximum number 
  of Search Events", y = "N", x = "Maximum number of 
       Search Events") + 
  theme_ipsum() + 
  theme(plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position="none", legend.box = "vertical")

bc# ggsave("transectsrepetitionsbx_p0_3-6.jpeg", dpi = 300)


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


# vl <- ggplot(data = sigma, aes(x=variable, y=value, fill = as.character(scenario))) 
# 
# vl + geom_violin(alpha = 1.2, scale = "width") +
#   geom_point(data = sigma_52, size = 3, alpha = 1.2,color = def_colors) +
#   facet_wrap(vars(stat), scales = "free_y",labeller = facet_lab) + 
#   scale_x_discrete(labels= x_labs) +
#   scale_fill_manual(values=my_colors) +
#   labs(title="σ - scale parameter", fill ="Search events", y = "sigma", x = "Number of Search Events") +
#   theme_ipsum() + 
#   theme(plot.title = element_text(size = 11),
#         axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
#         axis.text.y = element_text(size = 11),
#         legend.position="none", legend.box = "horizontal")
# 
# ggsave("transectsrepetitionsvl_sigma_3-6.jpeg", dpi = 300)



bx <- ggplot(data = sigma, aes(x=variable, y=value, fill = as.character(scenario))) 

bx + geom_boxplot(width = 0.6, alpha = 1.5) +
  geom_boxplot(data = sigma_52, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y",labeller=facet_lab) + 
  scale_x_discrete(labels= x_labs_t) +
  scale_fill_manual(values=my_colors) +
  labs(title="σ - scale parameter", fill ="Maximum number 
  of Search Events", y = "N", x = "Maximum number of 
       Search Events") + 
  theme_ipsum() + 
  theme(plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position="bottom", legend.box = "horizontal")

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
