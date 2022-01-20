#' @title Plot trace plots and posterior density of parameters
#'
#' @description
#' \code{PlotJagsParams} Plot trace plots, posterior density of parameters, and a vertical red line for the true parameter value (if known)
#'
#' @param jags.samples A \code{mcmc.list} object from coda.samples()functions in rjags
#' @param params A \code{vector} of strings with the variable names for plots should be computed
#' @param sim.values A \code{vector} with the true values of the parameters in \code{params}
#' @param trace.plot A \code{logial} if traceplot should be plotted
#' @param density.plot A \code{logial} if posterior density should be plotted
#'
#'
#' @return Plot trace plots and posterior density of parameters  
#' @example 
#' 




PlotJagsParams <- function( jags.samples 
                          , params = NULL
                          , sim.values = NULL
                          , trace.plot = TRUE
                          , density.plot = TRUE)
{
if(is.null(params)){params <- colnames(jags.samples[[1]])}
for(i in 1:length(params))
   {
   if(trace.plot & density.plot){par(mfrow=c(1,2))}
   if(trace.plot){
      traceplot(jags.samples[ ,params[i]])
      if(!is.null(sim.values)){abline(h = sim.values[i], col = "black", lwd = 3, lty = 2)}}
   if(density.plot){
      plot(density(unlist(jags.samples[ ,params[i]])), main = params[i])
      if(!is.null(sim.values)){abline(v = sim.values[i], col = "red", lwd = 2)}}
   }
}