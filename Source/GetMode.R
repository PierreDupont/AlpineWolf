

GetMode<-function(x,int=TRUE,plot.check=FALSE){
   
   if(!int){
   dens<-density(x,na.rm=TRUE)
   this.mode<-dens$x[dens$y==max(dens$y)][1]
   }else{
      this.mode=as.numeric(names(table(x))[table(x)==max(table(x))])[1]
     }
   if(plot.check){
   plot(dens)
   abline(v=this.mode,col="red")
   }
   return(this.mode)
}
   
x<-round(rnorm(1000,13,10),0)
# GetMode(x,plot.check=TRUE)   
# GetMode(rpois(1000,6),plot.check=TRUE)   

