#' @title Function to return the most common value in a list.
#'
#' @description
#' \code{GetMostCommon} returns a list object with \code{} 
#' 
#' @param data A vector.
#'
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/GetMostCommon.R
#' @keywords simul
#'
#' @examples
#'
#' GetMostCommon(x=c(1,2,3,4,NA,NA,NA,2,3))
#' 





GetMostCommon<-function(x,priority=NULL)
{
   if(all(is.na(x)))
   {
      temp<-NA
   }else
   {
      temp<-data.frame(table(x))
      temp<-temp[temp[,2]==max(temp[,2]),1][1]
      if(!is.null(priority))
      {	
         if(any(x%in%priority))
         {
            temp<-most(x[x%in%priority])
         }
      }
   }
   out<-as.character(temp)
   if(is.numeric(x))out<-as.numeric(out)
   return(out)
}