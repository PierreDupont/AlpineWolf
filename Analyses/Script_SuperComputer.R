## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())

## ------ IMPORT REQUIRED LIBRARIES ------
library(nimble)
library(nimbleSCR)

## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)


##---- Set-up directories
inPath <-"./input/"       ## Directory with input files as .RData
procPath <-"./processed/" ## Directory to store processed input files
outPath <-"./output/"     ## Directory to store processed output files


##---- Load one input file at random
file.list <- list.files(inPath)
set <- sample(file.list, 1)
load(paste(path, set, sep = ""))
print(set)

##---- Move the selected input file to the 'procPath' directory (to avoid repeated runs).
file.rename( from = file.path(inPath, set),
             to = file.path(procPath, set))                                 

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
                       thin = 1)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = nimModel, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc

##---- Run nimble MCMC in multiple bites
outFileName <- paste0()
print(system.time(
  runMCMCbites( mcmc = Cmcmc,
                bite.size = 500,
                bite.number = 40,
                path = file.path(outPath,set))))