## ------   4. FIT MODEL -----
## ------ SET REQUIRED WORKING DIRECTORIES ------

project_dir = "/archive/home/mboiani/AlpineWolf"
modelName = "AlpineWolf.5.4_RJMCMC_2324"
SourcePath <- "/archive/home/mboiani/AlpineWolf/Source"
## ------ IMPORT REQUIRED LIBRARIES ------
library(nimble)
library(nimbleSCR)
library(R.utils)

## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(SourcePath, modifiedOnly = F)


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
c_value <- c(1)
##---- Set-up directories
inPath <- file.path(modelName, "input")       ## Directory with input files as .RData
outPath <- file.path(modelName, "output")     ## Directory to store processed output files


for(c in c_value){
  file_to_load <- file.path(inPath, paste0(modelName,"_", c, ".RData"))
  if(!file.exists(file_to_load)) stop("File does not exist: ", file_to_load)
  load(file_to_load)
  
  ## Create nimble model object etc. here if it depends on the loaded data
  nimModel <- nimbleModel( code = modelCode,
                           constants = nimConstants,
                           inits = nimInits,
                           data = nimData,
                           check = FALSE,
                           calculate = FALSE) 
  nimModel$calculate()
  
  CsimModel <- compileNimble(nimModel)
  CsimModel$calculate()
  
  nimParams <- c("N", "p0", "sigma", "psi", "zRJ","psiRJ",
                 "betaDet", "betaHab.raw", "betaHab", "theta", "rho",
                 "z", "s", "status", "sex")
  
  conf <- configureMCMC( model = nimModel,
                         monitors = nimParams,
                         thin = 1)
  
  configureRJ( conf =  conf,
               targetNodes = 'betaHab.raw',
               indicatorNodes = 'zRJ',
               control = list(mean = 0, scale = .2))
  
  Rmcmc <- buildMCMC(conf)
  compiledList <- compileNimble(list(model = nimModel, mcmc = Rmcmc))
  Cmodel <- compiledList$model
  Cmcmc <- compiledList$mcmc
  
  print(system.time(
    runMCMCbites( mcmc = Cmcmc,
                  model = Cmodel,
                  bite.size = 1000,
                  bite.number = 100,
                  path = file.path(outPath, paste0(modelName, "_", c, ".RData")),
                  save.rds = TRUE)))
}