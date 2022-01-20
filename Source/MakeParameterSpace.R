
#' @title Function to Generate a dataframe of parameter space for simulation purposes 
#'
#' @description
#' \code{GenerateParameterSpace} returns a dataframe with the different parameters to be used in simulations 
#' 
#' @param list.param A \code{list} of vectors of the different parameters 
#' @param n.rep A \code{numeric}  Numeric variable denoting the number of repetition a specific combination should be repeated
#

#' @keywords simul
#'
#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' list.param <- list(  N = seq(5,30,by=5), sigma = c(5000, 12000, by=2000), p0 = c(1:5), beta = c(1))
#' GenerateParameterSpace(list.param, n.rep=10)   
#' 


MakeParameterSpace <- function(list.param = list.param
                                   ,n.rep = 10 
){
   param.df   <- expand.grid(list.param)
   n.comb <- nrow(param.df)
   param.df <- param.df[rep(seq_len(nrow(param.df)), n.rep), ]
   param.df$set_ID <- rep(seq_len(n.comb), n.rep)
   
   param.df <- param.df[sort(param.df$set_ID),]
   param.df$rep_ID <- rep(c(1:n.rep), n.comb)
   
   return(param.df)
}  
