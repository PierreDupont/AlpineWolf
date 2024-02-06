#' @title Run bite-size NIMBLE MCMC algorithms
#' 
#' @description
#' \code{runMCMCbites} is a wrapper R function to run a compiled NIMBLE model in 
#' multiple small size bites. This reduces memory usage and allows saving MCMC 
#' samples on the fly.
#' 
#' \code{restartMCMCbites} is a wrapper R function to restart and extend a previous
#' MCMC run, starting from the last point it was saved (as .rds file). In addition to
#' the saved model and MCMC states, this requires the MCMC configuration, nimble model
#' and MCMC to be built as in the original MCMC run.
#' By default, \code{restartMCMCbites} identifies and uses the last .rds file as a starting point.
#' the \code{path.rds} option allows to use a specific .rds file as starting point instead.
#' 
#' \code{collectMCMCbites} is a convenience R function to combine multiple bites
#' from the same MCMC chain into a single \code{mcmc} object (as defined in the 
#' \code{coda} package). If multiple chains were run, MCMC bites are combined 
#' chain-wise and compiled into a \code{mcmc.list} object.
#' 
#' @param mcmc a \code{NIMBLE MCMC algorithm}. See details. 
#' @param bite.size an \code{integer} denoting the number of MCMC iterations in 
#' each MCMC bite. 
#' @param bite.number an \code{integer} denoting the number of MCMC bites to be run.
#' @param path a \code{character string} of the path where MCMC bite outputs are
#'  to be saved as .RData files (when using \code{runMCMCbites}) or looked for 
#'  (when using \code{collectMCMCbites}). 
#' @param burnin an \code{integer} denoting the number of MCMC bites to be removed 
#' as burn-in. 
#' @param pattern a \code{character string} denoting the name of the object containing
#' MCMC samples in each .R Data bite file.
#' @param param.omit a \code{character vector} denoting the names of parameters 
#' that are to be ignored when combining MCMC bites.
#' @param progress.bar a \code{logical value}. if TRUE, a separate progress bar is 
#' printed for each MCMC chain.

##-- Run MCMC in multiple bites
runMCMCbites <- function( mcmc,
                          model = NULL,
                          bite.size,
                          bite.number,
                          path,
                          save.rds = FALSE){
  if(!dir.exists(path))dir.create(path, showWarnings = F, recursive = T)
  if(save.rds & is.null(model))stop("You must provide a nimble model object in order to use the 'save.rds' option")
  
  ##-- Loop over number of bites
  ptm <- proc.time()
  for(nb in 1:bite.number){
    print(nb)
    if(nb == 1){
      ##-- Run initial MCMC
      MCMCRuntime <- system.time(mcmc$run(bite.size))
    } else {      
      ##-- Run subsequent MCMCs
      MCMCRuntime <- system.time(mcmc$run(bite.size,
                                          reset = FALSE))
      ##-- Remove the last .rds file 
      if(save.rds)file.remove(file.path(path, paste0("MCMC_bite_", nb-1, ".rds")))   
    }
    
    ##-- Store bite output in a matrix
    mcmcSamples <- as.matrix(mcmc$mvSamples)
    mcmcSamples2 <- as.matrix(mcmc$mvSamples2)
    CumulRuntime <- proc.time() - ptm
    
    ##-- Export NIMBLE output 
    outname <- file.path( path, paste0("MCMC_bite_", nb, ".RData"))
    if(is.null(dim(mcmcSamples2))){
      save( CumulRuntime,
            MCMCRuntime,
            mcmcSamples,
            file = outname)
    } else {
      save( CumulRuntime,
            MCMCRuntime,
            mcmcSamples,
            mcmcSamples2,
            file = outname)
    }
    
    ##-- Save .rds file (current state of the model)
    if(save.rds){
      saveRDS(list( modelState = getModelState(model),
                    mcmcState = getMCMCstate(conf, mcmc)),
              file = file.path(path, paste0("MCMC_bite_", nb, ".rds")))
    }
    
    ##-- Free up memory space
    rm("mcmcSamples") 
    mcmc$mvSamples$resize(0)  
    mcmc$mvSamples2$resize(0) 
    gc()                      
  }#nb
}


##-- Restart MCMC run from saved model state
restartMCMCbites <- function( model,
                              conf,
                              mcmc,
                              bite.size,
                              bite.number,
                              path,
                              save.rds = TRUE,
                              path.rds = NULL){
  ##-- Load the .rds file 
  if(is.null(path.rds)){
    file.rds <- list.files(path)[grepl(".rds",list.files(path))]
    path.rds <- file.path(path, file.rds[length(file.rds)])
  }
  stateList <- readRDS(path.rds)
  
  ##-- Extract number of bites from the name of the .rds file 
  nbBites <- as.numeric(gsub(".*(_|\\s)(.*).rds", "\\2", path.rds)) 
  
  ##-- Restore the saved "state" into the new model and new MCMC
  setModelState(model, stateList$modelState)
  setMCMCstate(conf, mcmc, stateList$mcmcState)
  
  ##-- Loop over number of bites
  ptm <- proc.time()
  for(nb in (nbBites+1):(nbBites+bite.number)){
    print(nb)
    if(nb == 1){
      ##-- Run initial MCMC
      MCMCRuntime <- system.time(mcmc$run(bite.size))
    } else {      
      ##-- Run subsequent MCMCs
      MCMCRuntime <- system.time(mcmc$run(bite.size,
                                          reset = FALSE))
      ##-- Remove the last .rds file 
      if(save.rds)file.remove(file.path(path, paste0("MCMC_bite_", nb-1, ".rds")))   
    }
    
    ##-- Store bite output in a matrix
    mcmcSamples <- as.matrix(mcmc$mvSamples)
    mcmcSamples2 <- as.matrix(mcmc$mvSamples2)
    CumulRuntime <- proc.time() - ptm
    
    ##-- Export NIMBLE output 
    outname <- file.path( path, paste0("MCMC_bite_", nb, ".RData"))
    if(is.null(dim(mcmcSamples2))){
      save( CumulRuntime,
            MCMCRuntime,
            mcmcSamples,
            file = outname)
    } else {
      save( CumulRuntime,
            MCMCRuntime,
            mcmcSamples,
            mcmcSamples2,
            file = outname)
    }
    
    ##-- Save .rds file (current state of the model)
    if(save.rds){
      saveRDS(list( modelState = getModelState(model),
                    mcmcState = getMCMCstate(conf, mcmc)),
              file = file.path(path, paste0("MCMC_bite_", nb, ".rds")))
    }
    
    ##-- Free up memory space
    rm("mcmcSamples") 
    mcmc$mvSamples$resize(0)  
    mcmc$mvSamples2$resize(0) 
    gc()                      
  }#nb
}


##-- Combine MCMC bites 
collectMCMCbites <- function( path,
                              burnin = 0,
                              pattern = "mcmcSamples",
                              param.omit = NULL,
                              progress.bar = T){
  require(coda)
  
  ##-- Two possibilities:
  if(length(list.dirs(path, recursive = F)) == 0){
    ## 1 - the path contains multiple bite files (== one chain)
    path.list <- path
    outDir <- NULL
  } else {
    ## 2 - the path contains multiple directories(== multiple chains)
    ## List the directories containing bite outputs
    outDir <- list.files(path, pattern = ".RData", ignore.case = T)
    path.list <- file.path(path, outDir)
  }
  
  ##-- Retrieve the minimum number of MCMC bites per directory
  num.bites <- unlist(lapply(path.list, function(x){
    length(list.files(x, pattern = ".RData", ignore.case = T))
    }))
  num.bites <- min(num.bites)
  if(num.bites <= burnin)stop("Number of MCMC bites to burn is larger than the number of bites available")   
  
  ##-- Set-up progress bar 
  if(progress.bar){
    pb = txtProgressBar( min = burnin+1,
                         max = num.bites,
                         initial = 0,
                         style = 3) 
  }
  
  ##-- Loop over the different MCMC chains
  res <- res2 <- list()
  for(p in 1:length(path.list)){
    print(paste("Processing MCMC chain", p, "of", length(path.list)))
    
    ## List all MCMC bites in directory p
    out.files <- list.files(path.list[p], pattern = ".RData", ignore.case = T)
    
    ## Check the order of the MCMC bites
    newOrder <- order(as.numeric(gsub("[^\\d]+", "", out.files, perl=TRUE)))
    out.files <- out.files[newOrder]
    
    ## Loop over MCMC bites
    out.list <- out.list2 <- list()
    for(b in (burnin+1):num.bites){
      ## Load bite number "x"
      load(file.path(path.list[p],out.files[b]))
      
      ## Check that the mcmc samples object exists
      objInd <- which(ls() == pattern)
      if(length(objInd)<=0){stop(paste0("no object called ", pattern,
                                        " was found in ", outDir[p], "/", out.files[b],
                                        "!"))}
      ## Get the mcmc samples object
      out <- get(ls()[objInd])
      
      ## Remove parameters to ignore (optional) 
      paramSimple <- sapply(strsplit(colnames(out), split = '\\['), '[', 1)
      paramInd <- which(! paramSimple %in% param.omit)
      out.list[[b]] <- out[ ,paramInd] 
      
      ## Get the mcmc samples 2 object (only if thin2 used)
      objInd2 <- which(ls() == paste0(pattern,"2"))
      sp2yes <- length(objInd2 > 0)
      if(sp2yes){
        out2 <- get(ls()[objInd2])
        paramSimple2 <- sapply(strsplit(colnames(out2), split = '\\['), '[', 1)
        paramInd2 <- which(! paramSimple2 %in% param.omit)
        out.list2[[b]] <- out2[ ,paramInd2]   
      }
      ## Print progress bar
      if(progress.bar){ setTxtProgressBar(pb,b) }
    }#b
    if(progress.bar)close(pb)
    
    ## Combine into matrices
    out.mx <- do.call(rbind, out.list)
    res[[p]] <- as.mcmc(out.mx)
    if(sp2yes){
      out.mx2 <- do.call(rbind, out.list2)
      res2[[p]] <- as.mcmc(out.mx2)
    }
  }#p
  
  if(sp2yes){
    res <- as.mcmc.list(res)
    res2 <- as.mcmc.list(res2)
    return(list(samples = res,
                samples2 = res2))
  }else{
    res <- as.mcmc.list(res)
    return(res)
  }
  
}


##-- Get state variable names
getStateVariableNames <- function(samplerDef) {
  resetMethod <- body(samplerDef$reset)
  stateVars <- character()
  if(resetMethod[[1]] != '{') stop('something wrong')
  numLines <- length(resetMethod)
  for(i in 1:numLines) {
    if(i == 1) next
    thisLine <- resetMethod[[i]]
    if(thisLine[[1]] == '<<-') {
      LHS <- thisLine[[2]]
      if(!is.name(LHS)) stop('haven\'t dealt with non-name-LHS case yet')
      stateVars <- c(stateVars, as.character(LHS))
    }
    if('my_calcAdaptationFactor' %in% all.names(thisLine)) {
      stateVars <- c(stateVars, 'my_calcAdaptationFactor')
    }
  }
  setupMethod <- body(samplerDef$setup)
  if('empirSamp' %in% all.names(setupMethod)) stateVars <- c(stateVars, 'empirSamp')
  return(stateVars)
}


##-- Get & set model state
getModelState <- function( model) {
  modelVarNames <- model$getVarNames()
  modelVarValuesList <- vector('list', length(modelVarNames))
  names(modelVarValuesList) <- modelVarNames
  for(var in modelVarNames) {
    modelVarValuesList[[var]] <- model[[var]]
  }
  return(modelVarValuesList)
}

setModelState <- function( model,
                           modelState) {
  modelVarNames <- model$getVarNames()
  if(!identical(sort(modelVarNames), sort(names(modelState)))) stop('saved model variables don\'t agree')
  for(var in modelVarNames) {
    model[[var]] <- modelState[[var]]
  }
  invisible(model$calculate())
}

##-- Get & set MCMC state
getMCMCstate <- function( conf,
                          mcmc) {
  stateVarNamesList <- vector('list', length(conf$samplerConfs))
  mcmcStateValuesList <- vector('list', length(conf$samplerConfs))
  for(i in seq_along(conf$samplerConfs)) {
    samplerDef <- conf$getSamplerDefinition(i)
    theseStateNames <- getStateVariableNames(samplerDef)
    theseStateValuesList <- vector('list', length(theseStateNames))
    names(theseStateValuesList) <- theseStateNames
    for(j in seq_along(theseStateNames)) {
      if(is.nf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$timesAdapted,
                                            gamma1 = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$gamma1)
        } else
          theseStateValuesList[[j]] <- mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]
      }
      if(is.Cnf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted'),
                                            gamma1 = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1'))
        } else
          theseStateValuesList[[j]] <- valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], theseStateNames[j])
      }
    }
    mcmcStateValuesList[[i]] <- theseStateValuesList
  }
  return(mcmcStateValuesList)
}

setMCMCstate <- function( conf,
                          mcmc,
                          mcmcState) {
  if(length(mcmcState) != length(conf$samplerConfs)) stop('saved mcmc samplers don\'t agree')
  for(i in seq_along(conf$samplerConfs)) {
    theseStateValuesList <- mcmcState[[i]]
    for(j in seq_along(theseStateValuesList)) {
      samplerStateName <- names(theseStateValuesList)[j]
      if(is.nf(mcmc)) {
        if(samplerStateName == 'my_calcAdaptationFactor') {
          mcmc$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$timesAdapted <- theseStateValuesList[[j]]$timesAdapted
          mcmc$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$gamma1 <- theseStateValuesList[[j]]$gamma1
        } else {
          mcmc$samplerFunctions$contentsList[[i]][[samplerStateName]] <- theseStateValuesList[[samplerStateName]]
        }
      }
      if(is.Cnf(mcmc)) {
        if(samplerStateName == 'my_calcAdaptationFactor') {
          invisible(valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted', theseStateValuesList[[j]]$timesAdapted))
          invisible(valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1', theseStateValuesList[[j]]$gamma1))
        } else {
          invisible(valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], samplerStateName, theseStateValuesList[[samplerStateName]]))
        }
      }
    }
  }
}