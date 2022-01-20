#' @title Function to a Search Grid (actually search points) into a search trail that follows the topography. Requires the output from MakeSearchGrid.
#' #'
#' @description
#' \code{MakeTopoSearchTrail} returns a list object with \code{} 
#' 
#' @param detector.sp SpatialPointsDataframe with locations of search focal points to be visited. Object (data attribute) requires a field identifying main cell id (main.cell.id) that the subplot fall into. 
#' @param transitionLayer A transition layer object that defines transitions costs between raster cells (typically based on topography).
#' @param dem.r A raster of the digital elevation model on which  (\code{transitionLayer}) is based.
#' @param plot Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated execution.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/MakeTopoSearchTrail.R
#' @keywords simul
#'
#' @examples
#' # Make a nested search grid:
#' MakeTopoSearchTrail(...)
#' 
#' 
#'  


MakeTopoSearchTrail<-function(
   detector.sp,
   transitionLayer,
   dem.r=dem.area,
   plot=TRUE
){
   
   
   if(plot){
      #-----HILLSHADE
      slope <- terrain(dem.r, opt='slope')
      aspect <- terrain(dem.r, opt='aspect')
      hill <- hillShade(slope, aspect, 65, 270)
      plot(detector.sp,col="white",pch=19,cex=0.6)
      plot(hill, col=grey(0:100/100),add=TRUE)
      image(dem.r,col=terrain.colors(25, alpha=0.35), add=TRUE)
      plot(detector.sp,col="blue",pch=19,cex=0.6,add=TRUE)

   }
   
   
   main.cell.id.list<-sort(unique(detector.sp@data$main.cell.id))
   
   trail.list<-  lapply(main.cell.id.list,function(i){
      this.xy<-detector.sp[(detector.sp@data$main.cell.id==i),]
      #plot(this.xy,add=TRUE,col="orange",pch=19)
      this.xy<-coordinates(this.xy)
      
      if(length(this.xy[])>=4){
         
         #   RANDOM ORDER OF LINE SEGMENTS - could change to a fixed sequence, e.g. one that results in a circular route.
         this.xy<-this.xy[sample(length(this.xy[,1])),]
         
         AtoB<-list()
         for(j in 2:length(this.xy[,1])){
            print(j)
            
            outLC<-MakeLCPath(x=transitionLayer,this.xy[j-1,],this.xy[j,],plot=FALSE)
            plot(outLC$AtoB.sl,add=TRUE,col="blue")
            AtoB[[j-1]]<-outLC$AtoB.sl
            
         }
         out<-AtoB[[j-1]]
         if(length(AtoB)>1) out<-do.call(bind,AtoB)
         return(out)
      }
      
   })
   
   
   trail.list<-trail.list[!unlist(lapply(trail.list,is.null))]
   
   trail.sl <- do.call(bind, trail.list)  
   

   
   return(trail.sl)
   
   
   
}