#' @title Converting the scale parameter of the half-normal in SCR to home range size (circular).
#' #'
#' @description
#' \code{Sigma2Area} returns a list with ....
#' 
#' @param sigma Numeric variable indicating the scale parameter of the half-normal.
#' @param p Numeric variable, indicating the probability of home range utilization.
#' @param plot.check Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated execution.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/Sigma2Area.R
#' @keywords SCR data prep
#'
#' @examples
#' #Sigma2Area:
#' 
#'
#' 
#' 
#'  

Sigma2Area<-function(sigma,p=0.95){
   q<-qchisq(0.95,2)#--FOLLOWING ROYLE BOOK PAGE 136
   radius<-sigma*sqrt(q)
   area<-pi*radius^2
   out<-list(area=area,radius=radius,sigma=sigma)
   return(out)
}