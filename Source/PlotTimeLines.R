##------TIME LINE PLOT-----##
PlotTimeLines <- function( myData   
                         , myDead = NULL
                         , years = NULL
                         , TITLE = NULL
                         , cols = NULL
                         , myData2 = NULL
                         , myDead2 = NULL
                         , cols2 = NULL)
{ 
   if(is.null(years))
   {
   years <- c(min(na.omit(myData$Year)),max(na.omit(myData$Year)))
   if(!is.null(myData2)){years <- c(min(na.omit(c(myData$Year,myData2$Year))),max(na.omit(c(myData$Year,myData2$Year))))}
   }
   ymax <- round(length(unique(myData$Id)),-2)
   if(!is.null(myData2)){ymax <- round(length(unique(myData$Id))+length(unique(myData2$Id)),-2)}
   if(is.null(cols)){cols <- c("gray80",rgb(239/250,43/250,45/250),rgb(0/250,40/250,104/250))}
   first <- aggregate(Year ~ Id, data = myData, min)
   names(first)<- c("Id","year.first")
   last <- aggregate(Year ~ Id, data = myData, max)
   names(last)<- c("Id","year.last")
   first <- merge(first,last,by="Id",all=TRUE)
   first <- first[order(first$year.first,-(first$year.last)), ]
   
   plot(1~1 ,type="n",ylim=c(0,ymax),xlim=c(years[1],years[2]),axes=FALSE,xlab=" ",ylab="Number of Individuals")
   abline(v=(years[1]+0.5):(years[2]+0.5),lty=2,col="grey")
   axis(1,at=(years[1]:years[2]),labels=(years[1]:years[2]),tck=0)
   axis(2,at=seq(0,ymax,ymax/10),labels=seq(0,ymax,ymax/10))
   if(!is.null(TITLE)){title(TITLE)}
   
   tick <- 0
   id <- as.character(first$Id)
   for ( i in 1:length(id))
      {
      ## draw segments from first to last detections
      detect <- c(first$year.first[first$Id==id[i]],first$year.last[first$Id==id[i]])
      segments(detect[1]-0.5, tick, detect[2]+0.5, tick, col=cols[1],lwd=2)
      ## add a colour for each "capture"
      temp <- myData[myData$Id==id[i],]
      for(j in unique(temp$Year))
         {
         k <- length(temp$Year[temp$Year==j])
         segments(j-0.5, tick, j+0.5, tick,col=cols[2],lwd=2)
         }
      tick <- tick+1
      }
   
   ## add a symbol for known mortality 
   if(!is.null(myDead))
      {
      tick<- 0
      for ( i in 1:length(id))
         {
         if(sum(myDead$Id==id[i])>=1)
            {
            d <- myDead$Death[myDead$Id==id[i]][1]
            points(d+0.5, tick, pch=19, cex=0.7, col=cols[3])
            }
         tick <- tick + 1 
         }
      }
   
   if(!is.null(myData2))
      {
      if(is.null(cols2)){cols2 <- c(rgb(249/250,205/250,48/250),rgb(22/250,101/250,161/250),"black")}
      first <- aggregate(Year ~ Id, data = myData2, min)
      names(first) <- c("Id","year.first")
      last <- aggregate(Year ~ Id, data = myData2, max)
      names(last) <- c("Id","year.last")
      first <- merge(first,last,by="Id",all=TRUE)
      first <- first[order(first$year.first,-(first$year.last)), ]
      
      tick <- length(unique(myData$Id))+1
      id <- as.character(first$Id)
      for ( i in 1:length(id))
         {
         ## Draw segments from first to last detections
         detect <- c(first$year.first[first$Id==id[i]],first$year.last[first$Id==id[i]])
         segments(detect[1]-0.5, tick, detect[2]+0.5, tick, col=cols2[1],lwd=2)
         ## Add a colour for each number of samples
         temp <- myData2[myData2$Id==id[i],]
         for(j in unique(temp$Year))
            {
            k <- length(temp$Year[temp$Year==j])
            segments(j-0.5, tick, j+0.5, tick,col=cols2[2],lwd=2)
            }
         tick <- tick+1
         }
      
      ## Add a symbol for known mortality 
      if(!is.null(myDead2))
         {
         tick<- length(unique(myData$Id))+1
         for ( i in 1:length(id))
            {
            if(sum(myDead2$Id==id[i])>=1)
               {
               d <- myDead2$Death[myDead2$Id==id[i]][1]
               points(d+0.5, tick, pch=19, cex=0.7, col=cols2[3])
               }
            tick <- tick + 1 
            }
         }  
      }
}