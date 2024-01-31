

# ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())

# read the c parameter from command line
# Usage:
# Rscript run_model_0.3_RJMCM_x.R C_VALUE

args = commandArgs(trailingOnly = TRUE)
c_value = args[1]



## ------ IMPORT REQUIRED LIBRARIES ------
library(nimble)
library(nimbleSCR)
library(R.utils)

project_dir = "/archive/home/ofriard/projects/virginia-boiani/SC_model"

## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = "Source", modifiedOnly = F)

modelName = "AlpineWolf_SC.0.4"

##---- Set-up directories
inPath <- file.path(modelName, "input/")       ## Directory with input files as .RData
procPath <- file.path(modelName, "processed/") ## Directory to store processed input files
outPath <- file.path(modelName, "output/")     ## Directory to store processed output files


## ------ III. FIT NIMBLE MODEL -----


load( file.path(inPath, paste0(modelName, "_", c_value, ".RData")))

##---- Create the nimble model object
nimModel <- nimbleModel( code = modelCode,
                         constants = nimConstants,
                         inits = nimInits,
                         data = nimData,
                         check = FALSE,
                         calculate = FALSE)
nimModel$calculate()

##---- Compile the nimble model object to C++
CsimModel <- compileNimble(nimModel)
CsimModel$calculate()

##---- Configure and compile the MCMC object
conf <- configureMCMC( model = nimModel,
                       monitors = nimParams,
                       thin = 10)


Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = nimModel, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc

##---- Run nimble MCMC in multiple bites
outFileName <- paste0()
mcmcRuntime <- system.time(
  runMCMCbites( mcmc = Cmcmc,
                model=nimModel,
                bite.size = 1000,
                bite.number = 100,
                path = file.path(outPath, paste0(modelName, "_", c_value, ".RData")),
                save.rds = TRUE))

print(mcmcRuntime)


#}#c
