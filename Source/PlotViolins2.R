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


plot.violins2<-function(dat.list
                        , x
                        , at = NULL
                        , violin.width = 0.20
                        , invlogit = FALSE
                        , ylim = NULL
                        , col = "darkblue"
                        , cex = 0.5
                        , add = FALSE
                        , plot.ci = 1
                        , border.col = "black"
                        , alpha = 0
                        , fromto = NULL
                        , median = TRUE
                        , scale.width = FALSE
                        , lwd.violins = 1 
                        , lwd.CI = 1){
   


   
   if(is.null(at))at<-1:length(x)
   if(!add){
      if(invlogit & is.null(ylim))ylim<-c(0,1)#inv.logit(range(unlist(dat.list)))
      if(!invlogit & is.null(ylim))ylim<-range(unlist(dat.list))
      plot(1,1,type="n",ylim=ylim,xlim=c(min(at)-0.5,max(at)+0.5),axes=FALSE,xlab="",ylab="")
      axis(1,at=at)#, labels=x,lwd=0)
      axis(2)
   }
   
   i<-1
   
   #---GET THE VIOLIN-SPECIFIC SCALE IF SO INDICATED
   # amp<-unlist(lapply(1:length(x),function(i){
   #    max(density(dat.list[[i]])$y)
   # }))
   
   for(i in 1: length(x)){
      
      
      
      
      temp<-density(dat.list[[i]])
      if(!is.null(fromto)) temp<-density(dat.list[[i]],from=fromto[1],to=fromto[2])
      if(invlogit){
         temp$x<-inv.logit(temp$x)
         dat.list[[i]] <- inv.logit( dat.list[[i]])
      }
      
   
         
      #violin.width<-0.20 #in units x (group variable)
      
      #if(scale.width) scal.y<-DoScale(c(0,temp$y),0,violin.width[1]*amp[i]/max(amp))[-1]#--scaled to a portion of the year for plotting
      
      if(scale.width) scal.y<-DoScale(c(0,temp$y),0,violin.width*max(temp$y))[-1]
      if(!scale.width) scal.y<-DoScale(c(0,temp$y),0,violin.width)[-1]#--scaled to a portion of the year for plotting
      
      #if(!scale.width & length(violin.width)==length(x)) scal.y<-DoScale(c(0,temp$y),0,violin.width[i])[-1]#--scaled to a portion of the year for plotting
      #if(!scale.width & length(violin.width)==1) scal.y<-DoScale(c(0,temp$y),0,violin.width)[-1]#--scaled to a portion of the year for plotting
      
      poly.x<-c(at[i]-scal.y,(at[i]+scal.y)[length(scal.y):1])
      poly.y<-c(temp$x,temp$x[length(temp$x):1])#---the number
      
      
      
      #points(poly.y~poly.x,type="l")
      
      yv<-poly.y
      xv<-poly.x
      

      rgb.col<-as.vector(col2rgb(col)/255)
      polygon.col<-adjustcolor(col,alpha=alpha)#rgb(rgb.col[1],rgb.col[2],rgb.col[3],alpha)
      polygon(xv,yv,col=NA,border=border.col, lwd = lwd.violins)
     
      
      #if(!is.null(plot.ci)){
      lci<-quantile(dat.list[[i]],(1-plot.ci)/2)
      uci<-quantile(dat.list[[i]],1-(1-plot.ci)/2)
      
      yvci<-yv[yv>=lci & yv<=uci]
      xvci<-xv[yv>=lci & yv<=uci]
      polygon(xvci,yvci,col=polygon.col,border=border.col, lwd = lwd.CI)#polygon.col
     # }
      
      # if(median==TRUE){
      #   this.median<-ifelse(invlogit,inv.logit(median(dat.list[[i]])),median(dat.list[[i]]))
      # }else{
      #   this.median<-ifelse(invlogit,inv.logit(mean(dat.list[[i]])),mean(dat.list[[i]]))
      # }
      
      if(median==TRUE){
         this.median<-ifelse(invlogit,(median(dat.list[[i]])),median(dat.list[[i]]))
      }else{
         this.median<-ifelse(invlogit,(mean(dat.list[[i]])),mean(dat.list[[i]]))
      }
      points(this.median~at[i],pch=19,col="white",cex=cex)
      
   }#####
   
}


