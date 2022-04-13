#' @title Function to re-scale a numeric variable or vector (used inside UTM2Grid function)
#'
#' @description
#' \code{DoScale} returns a list object with \code{} 
#' 
#' @param data Numeric variable or vector to be re-scaled 
#' @param l Numeric variable denoting the lower bound of the new scale.
#' @param u Numeric variable denoting the upper bound of the new scale.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/DoScale.R
#' @keywords simul
#'
#' @examples
#' # Scale some numbers to a range between 0 and 1:
#' DoScale(data=c(1,2,3,4,NA,NA,NA,2,3),l=0,u=1)
#' 



DoScale<-function( data,
                   l = 0,
                   u = 1){
   x <- data + min(data, na.rm = TRUE)
    
   if(all(is.na(x))){
      scaled.x<-rep(NA,length(data))
   } else {
      min.x <- min(x,na.rm=TRUE)
      max.x <- max(x,na.rm=TRUE)
      scaled.x <- (u-l)*(x-min.x)/(max.x-min.x) + l
   }
   return(scaled.x)
}