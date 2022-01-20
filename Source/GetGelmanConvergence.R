#' @title Computes the gelmanConvergence values
#'
#' @description
#' \code{GetGelmanConvergence} Computes the gelmanConvergence values from the gelman.diag function from the coda package. It uses a moving window over the iteration to determine (from which iteration)
#' convergence has been reached. it only does from a mcmc object (coda.samples function in rjags) #' 
#'
#' @param jagsouput A \code{mcmc.list} object from coda.samples()functions in rjags
#' @param variable.names A \code{vector} of strings with the variable names for which gelman values should be computed 
#' @param moving.window A \code{numeric} for the length of the moving window
#' @param n.bins A \code{numeric} if the moving window should go sequentially through each iteration (n.bins=1) of it can go every nbins iterations. (use to speed up the function)
#' @param accept.value A \code{numeric} with the treshold of acceptance (<= acceptvalue)
#' @param confidence  transformobject, autoburnin, multivariate are arguments of the  gelman.diag() function. 
#' @param plot.check A \code{logical} whether gelman.plots should be computed 


#' @return A \code{list} object with convergence values 
#' $overall.gelman: output from gelman.diag function
#' $iteration.convergence: the iteration (representing the x value on the left sise of the moving window) for which convergence was reached (if all values within the moving window are <= accept.value).  
#' @example 
#' 
#'  Gelman <- GetGelmanConvergence( jagsouput=jagsouput
#'                                 , variable.names = params
#'                                 , accept.value = 1.01
#'                                 , moving.window = 500 )




GetGelmanConvergence <- function(   jagsouput = jagsouput1
                                  , variable.names = variable.names
                                  , moving.window = 20
                                  , n.bins = 1
                                  , accept.value = 1.01
                                  , confidence = 0.95
                                  , transform = FALSE
                                  , autoburnin = TRUE
                                  , multivariate = TRUE
                                  , plot.check = TRUE){
   
   
   ########## do the subset based on the parameter that we want to follow ##
   colnamesx <- colnames(jagsouput[[1]])
   param.pos <- NA
   
   for(i in 1:length(variable.names)){
      pos <- grep(variable.names[i],colnamesx)
      param.pos[(length(param.pos)+1):(length(param.pos)+length(pos))] <- pos
   }
   param.pos <- param.pos[!is.na(param.pos)]
   
   
   for( i in 1:length(jagsouput)){
      jagsouput[[i]] <- jagsouput[[i]][,param.pos]
   }
   
   ### now proceed the window of acceptance### 
   startx <- start(jagsouput)
   
   n.iterations <- end(jagsouput)-  moving.window
   
   
   
  
 

  #<- floor(((niter(jagsouput)/thin(jagsouput)))/n.bins) 
   
   last.iter <- c(seq(from = start(jagsouput), to = n.iterations, by = (n.bins * thin(jagsouput))))
   arr <- array(0, dim=c(length(param.pos),2 ,length(last.iter)  ))
   
   for(i in 1: length(last.iter) ){
      arr[,,i] <- round(gelman.diag(  window(jagsouput, start=last.iter[i], end = last.iter[i] + moving.window )
                                    , confidence = confidence
                                    , transform = transform
                                    , autoburnin = autoburnin
                                    , multivariate = FALSE)$psrf, digits = 3)
   }
   
   iteration.convergence <- apply(arr,1, function(x) suppressWarnings(min(which(x[1,] <= accept.value)))  )
   #iteration.convergence <- apply(arr,1, function(x) min(which(x[1,] <= accept.value))) 
   
   iteration.convergence[iteration.convergence==Inf] <- NA
   iteration.convergence <- data.frame(iteration.convergence)
   
   row.names(iteration.convergence ) <- colnames(jagsouput[[1]])
   
   if(plot.check==T){gelman.plot(jagsouput)    }
   
   
    overallgelman <- gelman.diag(jagsouput 
                              , confidence = confidence
                              , transform = transform
                              , autoburnin = autoburnin
                              , multivariate = FALSE)
   
    return( list( overall.gelman = overallgelman
           ,iteration.convergence = iteration.convergence
           ))
   
   
}




