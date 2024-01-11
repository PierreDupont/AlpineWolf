#' Truncated binomial distribution
#'
#' The \code{dtbinom} distribution is a truncated version of the binomial distribution.
#'
#' @name dtbinom
#'
#' @param x a scalar denoting the number of successes.
#' @param size a scalar with the number of trials (zero or more).
#' @param prob a scalar denoting the success probability. 
#' @param a and @param b, two scalar parameters denoting the left and right truncation boundaries (included).
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' @param n Number of observations. Only n = 1 is supported.
#'
#' @return The log-likelihood value associated with the truncated binomial observation.
#'
#' @author Pierre Dupont
#'
#' @import nimble
#' @importFrom stats dbinom rbinom
#'
#' @examples
#' ## define model code
#' code <- nimbleCode({
#'     
#'     v[1:(G+1)] ~ dmulti( size = K,
#'                      prob = p[]
#'     )
#'     p ~ dunif(0,1)
#'     p_vector[1:J] <- p
#'     y[1:J] ~ dbinom_vector(size = trials[1:J],
#'                            prob = p_vector[1:J])
#' })
#'  
#' ## simulate truncated binomial data
#' J <- 1000
#' trials <- sample(x = 10, size = J, replace = TRUE)
#' y <- rtbinom(J, size = trials, prob = 0.21)
#'  
#' constants <- list(J = J, trials = trials)
#'  
#' data <- list(y = y)
#'  
#' inits <- list(p = 0.5)
#'  
#' ## create NIMBLE model object
#' Rmodel <- nimbleModel(code, constants, data, inits)
#'  
#' ## use model object for MCMC, etc.
#'
#' @export
NULL

#' @rdname dtbinom
#' @export
dtbinom <- nimbleFunction(
  run = function( x = double(0),                   
                  size = double(0),              
                  prob = double(0), 
                  a = integer(0, default = 0),
                  b = integer(0, default = 0),
                  log = integer(0, default = 0)
  ){
    ##-- Define output type
    returnType(double(0))
    
    ##-- Define boundaries
    if(a <= 0) a <- 0 
    if(b <= 0) b <- size
    
    ##-- Two possibilities
    ##-- 1 - x < a or x > b ==> prob = 0
    ##-- 2 - a <= x <= b ==> dbinom(x,size,prob)/ sum(dbinom(a:b,size,prob)
    
    ##-- Case 1: x < a or x > b
    if(x < a | x > b){
      if(log) return(-Inf) else return(0)
    }
    
    ##-- Case 2: a <= x <= b 
    ##-- Compute P(a < X < b) for X ~ Binomial(size,prob)
    normCst <- 0
    
    for(i in a:b){
      thisProb <- dbinom( x = i,
                         size = size,
                         prob = prob,
                         log = 0)
      #sumProb <- sum(vecProb)
      normCst <- normCst + thisProb
    }#i
    logProb <- dbinom(x,size,prob,log = 1) - log(normCst)
    if(log) return(logProb) else return(exp(logProb))
  })



#' @rdname dtbinom
#' @export
rtbinom <- nimbleFunction(
  run = function( n = double(0, default = 1),                   
                  size = double(0),              
                  prob = double(0), 
                  a = integer(0, default = 0),
                  b = integer(0, default = 0),
                  log = integer(0, default = 0)
  ){
    ##-- Define output type
    returnType(double(0))

    ##-- Return type declaration
    if(n!=1){print("rtbinom only allows n = 1; using n = 1")}

    ##-- Define boundaries
    if(a <= 0) a <- 0 
    if(b <= 0) b <- size
    
    ##-- Case 2: a <= x <= b 
    ##-- Compute P(a < X < b) for X ~ Binomial(size,prob)
    for(i in a:b){
      thisProb <- dbinom( x = i,
                          size = size,
                          prob = prob,
                          log = 0)
      normCst <- normCst + thisProb
    }#i
    logProb <- dbinom(x,size,prob,log = 1) - log(normCst)
    
  return(y) 
  })



##-- Register distribution
registerDistributions(
  list(
    dtbinom = list(
      BUGSdist ='dtbinom(size, prob, a, b)',
      Rdist = c('dtbinom(size, prob, a, b)'),
      types = c('value = double(0)',
                'size = double(0)',
                'prob = double(0)',
                'a = double(0)',
                'b = double(0)'),
      discrete = TRUE,
      mixedSizes = TRUE,
      pqAvail = FALSE)
  ), verbose = F)

