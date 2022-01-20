#' @title Converting home range size (circular) to the scale parameter of the half-normal in SCR.
#' #'
#' @description
#' \code{Area2Sigma} returns a list with ....
#' 
#' @param area Numeric variable indicating the size of the home range. 
#' @param p Numeric variable, indicating the probability of home range utilization.
#' @param plot.check Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated execution.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/Area2Sigma.R
#' @keywords SCR data prep
#'
#' @examples
#' #Area2Sigma:
#' 
#
#' 
#' 
#'  

Area2Sigma<-function(area,p=0.95){
   q<-qchisq(p,2)#--FOLLOWING ROYLE BOOK PAGE 136
   radius<-sqrt(area/pi)
   sigma<-radius/sqrt(q)
   out<-list(area=area,radius=radius,sigma=sigma)
   return(out)
}