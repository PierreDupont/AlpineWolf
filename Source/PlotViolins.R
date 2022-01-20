#' @title Function to make a violin plot.
#'
#' @description
#' \code{plot.violins} creates a violin plot from a list of numerical vectors.
#' 
#' @param dat.list A list of numerical vectors. Each list item will be associated with its own violin, representing its distribution.
#' @param x A character or numerical vector indicating the labels associated with each item in dat.list.
#' @param at The location of each violin on the X axis.
#' @param add If the violins to be added to an existing plot (TRUE) or if a new plot is to be generated (FALSE)
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/plot.violins.R
#' @keywords simul
#'
#' @examples
#' # Generate a violin plot
#' plot.violins(dat.list=lapply(1:5,function(x)rnorm(10000,x,5)),x=1:5,at=NULL,invlogit=FALSE,ylim=NULL,col="darkblue",cex=0.5,add=FALSE)
#' 



plot.violins<-function(dat.list,x,at=NULL, violin.width=0.20,invlogit=FALSE,ylim=NULL,col="darkblue",cex=0.5,add=FALSE,plot.ci=NULL,border.col="black",alpha=0,fromto=NULL,median=TRUE){
   


   
   if(is.null(at))at<-1:length(x)
   if(!add){
      if(invlogit & is.null(ylim))ylim<-c(0,1)#inv.logit(range(unlist(dat.list)))
      if(!invlogit & is.null(ylim))ylim<-range(unlist(dat.list))
      plot(1~1,type="n",ylim=ylim,xlim=c(min(at)-0.5,max(at)+0.5),axes=FALSE,xlab="",ylab="")
      axis(1,at=at)#, labels=x,lwd=0)
      axis(2)
   }
   
   i<-2
   for(i in 1: length(x)){
      
      temp<-density(dat.list[[i]])
      if(!is.null(fromto)) temp<-density(dat.list[[i]],from=fromto[1],to=fromto[2])
      if(invlogit){
         temp$x<-inv.logit(temp$x)
      }
      
      #violin.width<-0.20 #in units x (group variable)
      
      scal.y<-DoScale(c(0,temp$y),0,violin.width)[-1]#--scaled to a portion of the year for plotting
      poly.x<-c(at[i]-scal.y,(at[i]+scal.y)[length(scal.y):1])
      poly.y<-c(temp$x,temp$x[length(temp$x):1])#---the number
      
      
      #points(poly.y~poly.x,type="l")
      
      yv<-poly.y
      xv<-poly.x
      

      rgb.col<-as.vector(col2rgb(col)/255)
      polygon.col<-rgb(rgb.col[1],rgb.col[2],rgb.col[3],alpha)
      polygon(xv,yv,col=polygon.col,border=border.col)
      if(median==TRUE){
      this.median<-ifelse(invlogit,inv.logit(median(dat.list[[i]])),median(dat.list[[i]]))
      }else{
      this.median<-ifelse(invlogit,inv.logit(mean(dat.list[[i]])),mean(dat.list[[i]]))
      }
      points(this.median~at[i],pch=19,col="white",cex=cex)
      
      if(!is.null(plot.ci)){
      lci<-quantile(dat.list[[i]],(1-0.95)/2)
      uci<-quantile(dat.list[[i]],1-(1-0.95)/2)
      lci.x<-temp$y[(temp$x-lci)==min(temp$x-lci)][1]
      uci.x<-temp$y[(temp$x-uci)==min(temp$x-uci)][1]
      lci.x<-DoScale(c(0,lci.x,temp$y),0,violin.width)[2]
      poly.lci.x<-c(at[i]-lci.x,at[i]+lci.x)
      segments(poly.lci.x[1],lci,poly.lci.x[2],lci,col="white",lwd=2)
      uci.x<-DoScale(c(0,uci.x,temp$y),0,violin.width)[2]
      poly.uci.x<-c(at[i]-uci.x,at[i]+uci.x)
      segments(poly.uci.x[1],uci,poly.uci.x[2],uci,col="white",lwd=2)
      }
   }#####
   
}


