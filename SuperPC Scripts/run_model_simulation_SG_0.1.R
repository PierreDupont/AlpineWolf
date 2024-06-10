
## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list = ls())

## ------ IMPORT REQUIRED LIBRARIES & FUNCTIONS ------
library(nimble)
library(R.utils)
## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = "../Source", modifiedOnly = F)
modelName = "simulation_SG.0.1"
thisDir <- file.path(simDir, "SpatialGroup", modelName)

##---- Set-up directories
inPath <- file.path("input")       ## Directory with input files as .RData
procPath <- file.path(modelName, "processed/") ## Directory to store processed input files
outPath <- file.path("output")     ## Directory to store processed output files


##---- Load one input file at random
file.list <- list.files(inPath)

if (length(file.list) == 0) {print("all models are processed"); exit()}

set <- sample(file.list, 1)
load(paste(inPath, set, sep = ""))
print(set)

##---- Move the selected input file to the 'procPath' directory (to avoid repeated runs).
file.rename( from = file.path(inPath, set),
             to = file.path(procPath, set)) 

## -----------------------------------------------------------------------------
## ------ III. MODEL FITTING ------
##---- Create the nimble model object
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits)
Cmodel <- compileNimble(model)
Cmodel$calculate()
modelConf <- configureMCMC( model, monitors = nimParams)
modelMCMC <- buildMCMC(modelConf)
CmodelMCMC <- compileNimble( modelMCMC, project = model)
system.time(nimOutput <- runMCMC( CmodelMCMC,
                                  niter = 20000,
                                  nchains = 3,
                                  samplesAsCodaMCMC = T))


  
## ------ SAVE OUTPUTS ------
save(nimOutput, file.path(outPath, paste0(tools::file_path_sans_ext(set),".RData")))
# pdf(file = file.path(outPath, paste0("SpatialCount_",r,".pdf" )),
#     width = 10, height = 8)
# plot(nimOutput)
dev.off()

