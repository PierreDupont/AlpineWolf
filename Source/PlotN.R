PlotN <- function( data = myPopSize.list$summary
                 , col = NULL
                 , obs.counts = NULL
                 , time.labels = NULL
                 , Nmax = NULL)
{
if(max(data$`0.975 %`)<=150)
   {
   y.max <- round(max(data$`0.975 %`),-1)+10
   } 
if(max(data$`0.975 %`)>150 & max(data$`0.975 %`)<=2500)
   {
   y.max <- round(max(data$`0.975 %`),-2)+100
   }
if(max(data$`0.975 %`)> 2500)
   {
   y.max <- round(max(data$`0.975 %`),-3)+1000
   }
   
if(!is.null(Nmax))
   {
   y.max <- Nmax
   }
   
axis2 <- seq(0,y.max,y.max/(y.max/10))

if (is.null(time.labels)){time.labels <- data$Year}
 
if(is.null(col)){col <- c(rgb(16/255,78/255,139/255,1),rgb(16/255,78/255,139/255,0.3))}    
       
plot(data$Year,data$Mean,type="n",pch=21,axes = FALSE,xlab="",ylab="",ylim=c(0,y.max))
axis(2, at = axis2, labels = axis2, cex.axis = 1.2, las=1, tck=0.01)
axis(1, at = data$Year, labels = time.labels, tck=0.01, cex.axis = 1.2, las = 1, hadj=0.5)
polygon(c(data$Year,rev(data$Year))
        , c(data$`0.025 %`,rev(data$`0.975 %`))
        , border=F, col=col[2])
points(data$Year, obs.counts, type="l", lty=2, lwd=1 , col="firebrick3")
points(data$Year, data$Mean, type="l", lwd=2, col=col[1])
mtext ("Years",1,2.2,cex=1.2,font=2)
mtext ("N",2,2.6,cex=1.2,font=2)
}
