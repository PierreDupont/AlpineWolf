#' @title Function to process jagsoutput from the simulations
#' @description
#' \code{ProcessJagsSimulation} returns a list object with statistics for all parameters. see Walther, B. A. and Moore, J. L. 2005.  Ecography  for info about formulas. it also returns the  number of ids detected and true number of alive ids
#' 
#' @param parm.df A dataframe coming the from the GenerateParameterSpace function with the list of parameters that were used for the simulation 
#' @param parameters A vector of string with the parameters for which mean/95CI should be returned 
#' @param name.file.in A string with the first part of the name of the file that was used to run the model
#' @param name.file.jags A string with the first part of the name of the JAGS output 
#' @param PathJags A string with the path of the directory containing the jagsoutout
#' @param PathIn A string with the path of the directory containing the data used to run the jagsmodel 
#'
#' @examples
#' 
#' JagsSimulationOutput <- ProcessJagsSimulation( parm.df =  parm.df1                      
#' , parameters = c("N", "sigma","beta.p")
#' , name.file.in = "SimDataSet"       
#' , name.file.jags = "JagsOutFORSimDataSet"
#' )

ProcessJagsSimulation <- function(  parameters.space =  parm.df1                      
                                    , parameters = c("N", "sigma","beta.p")
                                    , name.file.in = "SimDataSet"       
                                    , name.file.jags = "JagsOutFORSimDataSet"
                                    , PathJags = PathJags
                                    , PathIn = PathIn
                                    , resolution.raster = resolution.raster
                                    , alive.state = 1
                                    , type = "coda"
                                    , gelman.convergence=TRUE
                                    , N.not.in.buffer=FALSE
                                    , rename=NULL
                                    , renameto=NULL
                                    , moving.window = 500
                                    , n.bins = 1
                                    , accept.value = 1.1
                                    , other.parameters= NULL
                                    , params.omit = NULL
                                    , compute.statistics=TRUE
                                    , myHabitat.list=NULL
                                    ,psi.value.in.parm=FALSE
                                    , multiplefilesforchains=FALSE)
   
{
 
   parameters1 <- parameters
   
   # FOR EACH SIMULATIONS 
   for(i in 1:dim(parameters.space)[1]){
   # ==== I. LOAD THE FILE ====    
   parmss <- parameters.space[i,]
   print("##########################")
   print(i)
   
   # Set path file 
   this.file  <- paste(name.file.jags, parmss$set_ID, "Rep", parmss$rep_ID, ".RData" , sep = "")
   this.input <- paste(name.file.in, parmss$set_ID, "Rep", parmss$rep_ID, ".RData" , sep = "")
   
   if(multiplefilesforchains!=TRUE){
      if(length(grep(this.file, list.files(PathJags)))==1){
      load(paste(PathJags, this.file, sep = ""))
      load(paste(PathIn, this.input, sep=""))
      }else{
         next
      }
   
    }else{

       jagsoutput.list<- list()
      for(j in 1:3){ # WE RAN 3 CHAINS 
         #LOAD JAGSOUTFILES
         # CHECK IF FILE EXISTS 
         
         if(length(grep(paste( name.file.jags, parmss$set_ID, "Rep", parmss$rep_ID, "Clus",j,".RData" , sep = ""), list.files(PathJags)))==1){      
            load(paste(PathJags,"/", name.file.jags, parmss$set_ID, "Rep", parmss$rep_ID, "Clus",j,".RData" , sep = ""))
            load(paste(PathIn,"/", name.file.in, parmss$set_ID, "Rep", parmss$rep_ID, "Clus",j,".RData" , sep = ""))
            
            jagsoutput.list[[j]] <- jagsoutput[[1]]
               }else{next}
         
      }#j
      
        if(length(jagsoutput.list)==0){next}
   # IF FILES DOENST EXIST REMOVE FROM THE LIST 
   non_null_names <- which(!sapply(jagsoutput.list, is.null))
   if(length(non_null_names) > 0){
      jags.output <- jagsoutput.list[non_null_names]
      class(jags.output) <- "mcmc.list"
   }
   
   }
   
   
      # ---- 1. DEAL WITH DIFFERENT  JAGS PACKAGE ----     
   if(exists("jagsoutput", envir = environment())){ ## for people like me that doesnt put the right name...
      jags.output <- jagsoutput
   }
   if(exists("NimbleOutput", envir = environment())){ ## for people like me that doesnt put the right name...
      jags.output <- NimbleOutput
   }

   
   
   if(type == "coda"){
      ## replace param names 
      if(!is.null(rename)){
         for(r in 1:length(rename)){
            for(l in 1:length(jags.output)){
            colnames(jags.output[[l]]) <- gsub(rename[r], renameto[r], colnames(jags.output[[l]]))
            }
         }
      }
      
      
      Process <- ProcessCodaOutput( x = jags.output, DIC = FALSE
                                    , params.omit = params.omit, verbose = TRUE)
      sims.list <- Process$sims.list
      colnames.sims <- Process$colnames.sims
      
      if(gelman.convergence==TRUE){
      Convergence <- GetGelmanConvergence( jagsouput = jags.output
                            , variable.names = parameters1
                            , moving.window = moving.window
                            , n.bins = n.bins
                            , accept.value = accept.value
                            , plot.check = FALSE)$iteration.convergence
      }
   }
   
   if(type == "RunParallelRjags"){
      Process <- ProcessCodaOutput( x = jags.output$samples, DIC = FALSE
                                  , params.omit = params.omit, verbose = TRUE)
      sims.list <- Process$sims.list
      colnames.sims <- Process$colnames.sims
      if(gelman.convergence==TRUE){
      Convergence <- GetGelmanConvergence( jagsouput = jags.output$samples
                                           , variable.names = parameters1
                                           , moving.window = moving.window
                                           , n.bins = n.bins
                                           , accept.value = accept.value
                                           , plot.check = FALSE)$iteration.convergence
      }
   }
   if(type == "R2jags"){
      sims.list <- jags.output$BUGSoutput$sims.list
   }
   if(type == "jagsUI"){
      Process <- ProcessCodaOutput( x = jags.output$samples, DIC = FALSE
                                      , params.omit = params.omit, verbose = TRUE)$sims.list
      
      sims.list <- jags.output$sims.list
      colnames.sims <- names(sims.list)
      if(gelman.convergence==TRUE){
         Convergence <- GetGelmanConvergence( jagsouput = jags.output$samples
                                              , variable.names = parameters1
                                              , moving.window = moving.window
                                              , n.bins = n.bins
                                              , accept.value = accept.value
                                              , plot.check = FALSE)$iteration.convergence
      }
   }
   
# ==== II. STORE INFORMATION  ====  
   if(N.not.in.buffer==TRUE){
      
      parmss$N.not.in.buffer  <- length(mySimulatedACs[!is.na(over(mySimulatedACs, myHabitat.list$habitat.poly)),])
      
      ## if density maps should be plotted #
      # myPopDensity <- GetPopDensity( habitat.mx = myHabitat.list$habitat.mx
      #                                , r.origin = myHabitat.list$habitat.r
      #                                , posterior.sy = sims.list$sxy[,,2]
      #                                , posterior.sx = sims.list$sxy[,,1]
      #                                , posterior.z = sims.list$z
      #                                , alive.states = 1 )
      
      myPopDensity.cut.list <-  EstimateN (r.origin = myHabitat.list$habitat.r             # A raster file of the study area 
                                           , IDCells.mx = myHabitat.list$IDCells.mx           # A matrix with the cells ID of the study area  
                                           , posterior.sy =  sims.list$sxy[,,2]          # The posterior x coordinates of individual ACs 
                                           , posterior.sx = sims.list$sxy[,,1]            # The posterior y coordinates of individual ACs
                                           , posterior.z = sims.list$z       # The posterior of z 
                                           , alive.states = 1           # Alive states in z
                                           , poly = myHabitat.list$habitat.poly                  # Polygon of the study area
                                           , time.interval = NULL     # Years to retained 
                                           , CI = c(0.025,0.975))$N.counts[,1]
      
      
      sims.list[["N.not.in.buffer"]] <- myPopDensity.cut.list
      
      df <- data.frame(iteration.convergence=NA)
      row.names(df) <- "N.not.in.buffer"
      Convergence <- rbind(Convergence, df)
      
      Parameters <- c(parameters1, "N.not.in.buffer")
   }else{
      Parameters <- c(parameters1)
   }
      # ---- 1. FOR THE FIRST FILE/SIMULATION CREATE EMPTY LIST  ----     
   
   if(i==1){ # 
      param.list <- list() 
      param.list[[Parameters[1]]] <- list()
      
      stats <- c("mean", "median", "sd", "lci", "uci", "range", "coverage","iteration.convergence","nb.detected","true.z.nb.alive", other.parameters) 
      if(N.not.in.buffer==FALSE){
         stats <- c(stats, "N.not.in.buffer") 
      }   
      if(compute.statistics==TRUE){
         stats1 <- c("ME", "SME", "CV", "RMSE","SRMSE") 
         stats <- c(stats, stats1)
      }
      for(j in 1:length(Parameters)){
         colnames.sims.j <- colnames.sims[grep(Parameters[j],colnames.sims)]
         
         for(s in 1:length(stats)){
            if(stats[s]!="range"){
               param.list[[Parameters[j]]][[stats[s]]] <- data.frame(X=0)
               colnames(param.list[[Parameters[j]]][[stats[s]]]) <- colnames.sims.j
               
            }else{
               param.list[[Parameters[j]]][[stats[s]]] <- data.frame(X=0,Y=0) 
               colnames(param.list[[Parameters[j]]][[stats[s]]]) <- colnames.sims.j
               
            }
         }
      }
   }

   
   # SIMULATION NAME 
   this.sim.name <- paste("set.", parmss$set_ID, "Rep.", parmss$rep_ID, sep="")

   
  if(psi.value.in.parm!=TRUE){
   if(sum(Parameters %in% "psi")>0){
      parmss$psi <- true.psi
   }
  }
      

      
      
      
      # ---- 2. FOR STORE STATISTICS   ----     
   
   for(j in 1:length(Parameters)){
      # Mean
      param.list[[Parameters[j]]][["mean"]][i,] <-  data.frame(t(colMeans(data.frame(sims.list[[Parameters[j]]] ))))
      rownames( param.list[[Parameters[j]]][["mean"]] )[i] <- this.sim.name
      # Median
      param.list[[Parameters[j]]][["median"]][i,] <-  data.frame(t(apply(data.frame(sims.list[[Parameters[j]]] ), 2,function(x) median(x))))
      rownames( param.list[[Parameters[j]]][["median"]] )[i] <-this.sim.name
      # SD
      param.list[[Parameters[j]]][["sd"]][i,] <-  data.frame(t(apply(data.frame(sims.list[[Parameters[j]]] ),2,function(x) sd(x))))
      rownames( param.list[[Parameters[j]]][["sd"]] )[i] <- this.sim.name
      # LCI
      param.list[[Parameters[j]]][["lci"]][i,] <-  data.frame(t(data.frame(apply(data.frame(sims.list[[Parameters[j]]]) , 2, function(x) quantile(x, 0.025)))))
      rownames(param.list[[Parameters[j]]][["lci"]])[i] <- this.sim.name
      # UCI
      param.list[[Parameters[j]]][["uci"]][i,] <-  data.frame(t(data.frame(apply(data.frame(sims.list[[Parameters[j]]]) , 2, function(x) quantile(x, 0.975)))))
      rownames(param.list[[Parameters[j]]][["uci"]])[i] <- this.sim.name
      # RANGE
      param.list[[Parameters[j]]][["range"]][i,] <-  data.frame(t(data.frame(apply(data.frame(sims.list[[Parameters[j]]]) , 2, function(x) range(x)))))
      rownames(param.list[[Parameters[j]]][["range"]])[i] <- this.sim.name

      param.list[[Parameters[j]]][["coverage"]][i,]  <-  data.frame( param.list[[Parameters[j]]][["lci"]][i,] <  parmss[,Parameters[j]]  &
                                                                   param.list[[Parameters[j]]][["uci"]][i,] >  parmss[,Parameters[j]] )
      rownames(param.list[[Parameters[j]]][["coverage"]])[i] <- paste("set.", parmss$set_ID, "Rep.", parmss$rep_ID, sep="")

      
      # ---- 3. FOR STORE STATISTICS COMPARING ESTIMATED / TRUE VALUES   ----     
      
      if(compute.statistics==TRUE){
      # All statistics can be found in Walther, B. A. and Moore, J. L. 2005. ECOGRAPHY
      # ME MEAN ERROR (bias)
  
         
      param.list[[Parameters[j]]][["ME"]][i,] <-  data.frame(mean(sims.list[[Parameters[j]]]) - parmss[,Parameters[j]])
      rownames(param.list[[Parameters[j]]][["ME"]])[i] <- this.sim.name
      
      # SME SCALED MEAN ERROR (relative bias) 
      param.list[[Parameters[j]]][["SME"]][i,] <-  data.frame(param.list[[Parameters[j]]][["ME"]][i,] / parmss[,Parameters[j]])
      rownames(param.list[[Parameters[j]]][["SME"]])[i] <- this.sim.name
      
      ## CV
      param.list[[Parameters[j]]][["CV"]][i,] <-   data.frame((100 * param.list[[Parameters[j]]][["sd"]][i,]) / param.list[[Parameters[j]]][["mean"]][i,])
      rownames(param.list[[Parameters[j]]][["CV"]])[i] <- this.sim.name
      
      # RMSE
      param.list[[Parameters[j]]][["RMSE"]][i,] <-   data.frame(sqrt( mean( (sims.list[[Parameters[j]]] - parmss[,Parameters[j]])^2) ) )
      rownames(param.list[[Parameters[j]]][["RMSE"]])[i] <- this.sim.name
      
      # SRMSE SCALED RMSE
      param.list[[Parameters[j]]][["SRMSE"]][i,] <-   param.list[[Parameters[j]]][["RMSE"]][i,] / parmss[,Parameters[j]] 
      rownames(param.list[[Parameters[j]]][["SRMSE"]])[i] <- this.sim.name
      
      }
   if(gelman.convergence==T){
      if(type != "R2jags"){
       param.list[[Parameters[j]]][["iteration.convergence"]][i,] <- data.frame(Convergence[Parameters[j],1])
       rownames(param.list[[Parameters[j]]][["iteration.convergence"]])[i] <- this.sim.name
      }
   }
   }#j
   
   
   
   # ---- 4. FOR STORE OTHER Parameters CONCERNING INDIVIDUALS  ----     
   
   if(sum(Parameters %in% "true.z.nb.alive") > 0){
      param.list[["true.z.nb.alive" ]][,i] <- data.frame(t(data.frame(apply(true.z, 2, function(x) sum(x %in% alive.state)) )))  ## carefull here state alive =2
      rownames(param.list[["true.z.nb.alive"]])[i] <- this.sim.name
   }
   
   if(sum(Parameters %in% "nb.detected")>0){
      param.list[["nb.detected"]][,i] <- data.frame(t(data.frame(apply(data.frame(my.jags.input$z), 2, function(x) sum(x %in% alive.state)) )))## carefull here state alive =2
      rownames(param.list[["nb.detected"]])[i] <- this.sim.name
   }
   
   
   
   if(!is.null(other.parameters) ){
      for(o in 1:length(other.parameters)){
         if(i==1){ 
         param.list[[other.parameters[o]]] <- data.frame(X=0)
         colnames( param.list[[other.parameters[o]]] ) <- other.parameters[o]
         }
         
            if( other.parameters[o] %in% c("y.all", "y.origin")){
               param.list[[other.parameters[o]]][i,] <- data.frame(sum(get0(other.parameters[o], envir = environment())))
               rownames(param.list[[other.parameters[o]]])[i] <- paste("set.", parmss$set_ID, "Rep.", parmss$rep_ID, sep="")
            }else{
               if(other.parameters[o] %in% c("my.jags.input")){
                  param.list[[other.parameters[o]]][i,] <- data.frame(nrow(get0("my.jags.input", envir = environment())$y))
                  rownames(param.list[[other.parameters[o]]])[i] <- this.sim.name
               }else{
                  if(other.parameters[o] %in% c("Time")){
                     param.list[[other.parameters[o]]][i,] <- data.frame(get0("Time", envir = environment())[3])
                     rownames(param.list[[other.parameters[o]]])[i] <- this.sim.name
               }
            }
         }
      }
    
   }
}#i
   
   # CONVERT ALL SIGMA ESTIMATES WITH THE RESOLUTION 
   if(sum(Parameters %in% "sigma") > 0 ){
      param.list[["sigma"]][1:5] <-  lapply( param.list[["sigma"]][1:5], function(x)  x * resolution.raster )  
   }#if
   
   return(param.list)
   }
   


