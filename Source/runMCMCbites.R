#' @title Run bite-size NIMBLE MCMC algorithms
#' 
#' @description
#' \code{runMCMCbites} is a wrapper R function to run a compiled NIMBLE model in 
#' multiple small size bites. This reduces memory usage and allows saving MCMC 
#' samples on the fly.
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
#' @param param.omit a \code{logical value}. if TRUE, a separate progress bar is 
#' printed for each MCMC chain.

runMCMCbites <- function( mcmc,
                          bite.size,
                          bite.number,
                          path){
  if(!dir.exists(path))dir.create(path)
  ptm <- proc.time()
  ## Loop over number of bites
  for(nb in 1:bite.number){
    print(nb)
    if(nb == 1){
      ## run initial MCMC
      MCMCRuntime <- system.time(mcmc$run(bite.size))
    } else {      
      ## run subsequent MCMCs
      MCMCRuntime <- system.time(mcmc$run(bite.size,
                                           reset = FALSE))
    }
    
    ## STORE BITE OUTPUT IN A MATRIX
    mcmcSamples <- as.matrix(mcmc$mvSamples)
    mcmcSamples2 <- as.matrix(Cmcmc$mvSamples2)
    CumulRuntime <- proc.time() - ptm
    
    ## EXPORT NIMBLE OUTPUT 
    outname <- file.path( path,
                          paste0("MCMC_bite_", nb, ".RData"))
    if(is.null(dim(mcmcSamples2))){
      save( CumulRuntime,
            MCMCRuntime,
            mcmcSamples,
            file = outname)
    }else{
      save( CumulRuntime,
            MCMCRuntime,
            mcmcSamples,
            mcmcSamples2,
            file = outname)
    }

    
    ## FREE UP MEMORY SPACE 
    rm("mcmcSamples") 
    mcmc$mvSamples$resize(0) ## reduce the internal mvSamples object to 0 rows
    mcmc$mvSamples2$resize(0) ## reduce the internal mvSamples object to 0 rows,
    gc() ## run R's garbage collector
  }#nb
}

collectMCMCbites <- function( path,
                              burnin = 0,
                              pattern = "mcmcSamples",
                              param.omit = NULL,
                              progress.bar = T){
  require(coda)
  
  ## Two possibilities:
  if(length(list.dirs(path, recursive = F)) == 0){
    ## 1 - the path contains multiple bite files (== one chain)
    path.list <- path
    outDir <- NULL
  } else {
    ## 2 - the path contains multiple directories(== multiple chains)
    ## List the directories containing bite outputs
    outDir <- list.files(path, ignore.case = T)
    path.list <- file.path(path, outDir)
  }
  
  ## Retrieve the minimum number of MCMC bites per directory
  num.bites <- unlist(lapply(path.list, function(x)length(list.files(x))))
  num.bites <- min(num.bites)
  if(num.bites <= burnin)stop("Number of MCMC bites to burn is larger than the number of bites available")   
  
  ## Set-up progress bar 
  if(progress.bar){
    pb = txtProgressBar( min = burnin+1,
                         max = num.bites,
                         initial = 0,
                         style = 3) 
  }
  
  ## Loop over the different MCMC chains
  res <- res2 <- list()
  for(p in 1:length(path.list)){
    print(paste("Processing MCMC chain", p, "of", length(path.list)))
    
    ## List all MCMC bites in directory p
    out.files <- list.files(path.list[p])
    
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