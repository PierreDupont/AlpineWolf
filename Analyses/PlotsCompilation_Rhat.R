
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
library(ggh4x)



## ------ 3. SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ 4.  SOURCE CUSTOM FUNCTIONS ------
# sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))


## -----------------------------------------------------------------------------
## ------ 5. LOAD RESULTS -----
## MODEL NAME 
modelName = "AlpineWolf.SC.1.1.0"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
res00 <- res
N00 <- res00$Rhat


modelName = "AlpineWolf.SC.1.1.1"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
res01 <- res
N01 <- res01$Rhat


modelName = "AlpineWolf.SC.1.1.4"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
res02 <- res
N02 <- res02$Rhat


modelName = "AlpineWolf.SC.1.1.2_A"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
res03 <- res
N03 <- res03$Rhat

modelName = "AlpineWolf.SC.1.1.2_B"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
res06 <- res
N06 <- res06$Rhat

modelName = "AlpineWolf.SC.1.1.2"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
res07 <- res
N07 <- res07$Rhat

modelName = "AlpineWolf.SC.1.1.5"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
res04 <- res
N04 <- res04$Rhat


modelName = "AlpineWolf.5.2_RJMCMC"
thisDir <- file.path(analysisDir, modelName)
load(file.path(thisDir, paste0(modelName, "_mcmc.RData")))
res05 <- res
N05 <- res05$Rhat


# Combine your 6 Rhat lists into one named list
model_list <- list(
  M1 = N00,
  M2 = N01,
  M3 = N02,
  M4 = N03,
  M5 = N04,
  M4a = N06,
  M4b = N07
  )

# Mapping raw parameter names to nice symbols
param_labels <- c( D = "Density", 
                   N = "Abundance", 
                   betaHab = "β", 
                   lambda0 = "λ", 
                   psi = "ψ", 
                   sigma = "σ", 
                   rho = "ρ", 
                   theta = "θ" )


# Modified flattening function with cleaner names 

flatten_model <- function(model, model_name) { 
  out <- purrr::imap_dfr(model, function(val, param) { 
    if (is.vector(val) || is.matrix(val)) { 
      val <- as.vector(val) 
      tibble(parameter = paste0(param, "_", seq_along(val)), Rhat = val) 
    } else { 
        tibble(parameter = param, Rhat = as.numeric(val)) 
    } 
  }) 

  out$model <- model_name 
    # Clean parameter labels 
  
out <- out %>% mutate(param_base = gsub("_\\d+", "", parameter), 
                        parameter_clean = case_when( param_base %in% 
                                                       names(param_labels) ~ paste0(param_labels[param_base], 
                                                                                    " (", parameter, ")"), TRUE ~ parameter )) 
return(out) 
}



all_models_df <- imap_dfr(model_list, flatten_model)


# Extract the symbol before the bracket and only append an index if repeated
all_models_df <- all_models_df %>%
  mutate(parameter_label = gsub("\\s*\\(.*\\)", "", parameter_clean)) %>%
  group_by(model, parameter_label) %>%
  mutate(parameter_pretty = if (n() > 1) {
    paste0(parameter_label, "_", seq_along(parameter_label))
  } else {
    parameter_label
  }) %>%
  ungroup()


all_models_df <- all_models_df %>%
  filter(!grepl("^(zRJ|BetaHab\\.Raw|BetaDet)", parameter))


all_models_df <- all_models_df %>%
  mutate(flag = ifelse(Rhat > 1.10, "Rhat > 1.10", "Rhat ≤ 1.10"))

# Split the data
df_small <- all_models_df %>% filter(model != "M5")
df_big   <- all_models_df %>% filter(model == "M5")

# p1: models M1–M4, stacked for clarity
p1 <- ggplot(df_small, aes(x = parameter_pretty, y = Rhat, fill = flag)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~model, scales = "free_y", ncol = 2) +
  # scale_y_continuous(limits = c(0,1.5), breaks = seq(0, 1.5, by = 0.5)) +
  scale_fill_manual(values = c("Rhat > 1.10" = "cadetblue", "Rhat ≤ 1.10" = "grey70")) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = "Parameter", y = "Rhat", title = NULL)

# p2: model M5, wider
p2 <- ggplot(df_big, aes(x = parameter_pretty, y = Rhat, fill = flag)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~model, scales = "free_y") +
  scale_fill_manual(values = c("Rhat > 1.10" = "cadetblue", "Rhat ≤ 1.10" = "grey70")) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right") +
  labs(x = "", y = "Rhat", title = NULL)

# Combine with relative widths (1 for p1, 2 for p2)
p1 + p2 + plot_layout(widths = c(1, 2))

