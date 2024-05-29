
## -------------------------------------------------------------------------- ##
## ------------------------ ALPINE WOLF SC ---------------------------------- ##
## -------------------------------------------------------------------------- ##
## ------1. CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------ 2. LOAD LIBRARIES ------
library(dplyr)
library(stringr)
library(abind)
library(R.utils)
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
library(patchwork)
library(mefa)


## ------ 3. SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ 4.  SOURCE CUSTOM FUNCTIONS ------
# sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))


def_colors <- ghibli::ghibli_palette("PonyoMedium")[4]
my_colors <- wes_palette("GrandBudapest2")[1:2]
my_colors1 <- wes_palette("GrandBudapest2")[4]

my_colors2 <- ghibli::ghibli_palette("PonyoMedium")[3:4]

mycolors <- c(my_colors, my_colors1, my_colors2)

## -----------------------------------------------------------------------------
## ------ 5. LOAD RESULTS -----
## MODEL NAME 
modelName = "AlpineWolf.SC.1.1.0"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_posterior.RData")))
res00 <- WA_Italy
N00 <- res00$PosteriorAllRegions
N_00 <- data.frame(N = N00,
                   model = "A.W.0")

modelName = "AlpineWolf.SC.1.1.1"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_posterior.RData")))
res01 <- WA_Italy
N01 <- res01$PosteriorAllRegions
N_01 <- data.frame(N = N01,
                   model = "A.W.1")

modelName = "AlpineWolf.SC.1.1.2"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_posterior.RData")))
res02 <- WA_Italy
N02 <- res02$PosteriorAllRegions
N_02 <- data.frame(N = N02,
                   model = "A.W.2")

modelName = "AlpineWolf.SC.1.1.4"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_posterior.RData")))
res03 <- WA_Italy
N03 <- res03$PosteriorAllRegions
N_03 <- data.frame(N = N03,
                   model = "A.W.3")

modelName = "AlpineWolf.SC.1.1.5"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_posterior.RData")))
res04 <- WA_Italy
N04 <- res04$PosteriorAllRegions
N_04 <- data.frame(N = N04,
                   model = "A.W.4")


N_all <- rbind(N_00, N_01, N_02, N_03, N_04)
N_original <- 952
## ------ 6. PLOT N -----
res_bias_ntot_plot <- N_all %>%
  # mutate(scenario = fct_relevel(scenario,"25", "50","75","100")) %>%
  ggplot(aes(x=model, y=N,fill=model)) +
  geom_violin(trim = FALSE) +
  geom_hline(aes(yintercept= N_original), color = "plum3", size=1) +
  geom_hline(aes(yintercept=816), color = "steelblue3", size=1) +
  geom_hline(aes(yintercept=1120), color = "steelblue3", size=1) +
  stat_summary(fun = median, geom = "point", color = "white", size = 0.5,
               position = position_dodge(0.9)) +
  # ylim(-1,1.5) +
  scale_fill_manual(values= wes_palette("Zissou1")
                    # , name = "Spatial Subsampling",labels=c('25%', '50%', '75%', '100%')
                    ) +
  labs(y = "Population size", x = "") +
  theme(legend.position = "right",
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        axis.text.x = element_text(size=0),
        axis.text.y = element_text(size=10))


res_bias_ntot_plot
# +   
#   geom_ribbon(aes(ymin = N_original - 136, ymax = N_original + 168), fill = "steelblue3") +
#   geom_line(aes(y = N_original)) 

# 816 952 1120





