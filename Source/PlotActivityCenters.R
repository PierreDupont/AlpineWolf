
#' Plot activity centers and link them with segments
#'
#'
#'
#' @author Cyril Milleret
#' @param sp a SpatialPointsDataFrame object with multiple points (or not) for each individiduals #
#'
#' @param ID a vector with IDs. Typically, a column of the sp object with the ID of each individuals
#' @param col.detections colors of the detection points.  "random" if colors should be selected randomly. otherwise provide a vector with color for each IDS e.g. c("black,"red","green")
#' @param col.ACS colors of the Acitivity centers
#' @param col.segments colors of the segments
#' @param cex size of the symbols
#' @param plot.segments Whether the segments should be plotted (TRUE) or not (FALSE)
#' @param plot.ACS Whether the ACs should be plotted (TRUE) or not (FALSE)
#'
#' @return Plot the activity centers and the associated points linked by a segment.
#' @usage ActivityCenters(sp=puechabonsp$relocs, ID=puechabonsp$relocs$Name, col.detections="random" ,
#' col.ACS="random", col.segments="random", cex=cex, plot.segments=TRUE, plot.ACS=TRUE )
#' @export
#' @examples
#'
#'
#' library(adehabitatHS)
#' data("puechabonsp")
#' image(puechabonsp$map)
#' points(puechabonsp$relocs)
#' ActivityCenters(sp=puechabonsp$relocs, ID=puechabonsp$relocs$Name, col.detections="random" ,
#' col.ACS="random", col.segments="random", cex=1, plot.segments=TRUE, plot.ACS=TRUE )
#'
#'


PlotActivityCenters <- function(sp=sp, ID=ID, col.detections="random" , col.ACS="random", col.segments="random",
                            cex=1, plot.segments=TRUE, plot.ACS=TRUE, pch.detect=16,bg.detect="black",lwd.segm=1){

  ID1 <- unique(ID)
  col <-  sample(colours(), length(ID1),replace = T)
  ## colors if random, otherwise use your own set of colors from a vector
  if(col.detections=="random"){
    col.detections = col
  }else{col.detections=col.detections}

  if(col.ACS=="random"){
    col.ACS = col
  }else{col.ACS=col.ACS}

  if(col.segments=="random"){
    col.segments = col
  }else{col.segments=col.segments}

for (i in 1:length(ID1)){## for each IDS ##

  tmp <- sp[which(ID==ID1[i]),]
  if(length(tmp)==1){# if 1 sample= just a point for the detection process
    points(tmp, pch=pch.detect, col=col.detections[i], cex=cex, bg=bg.detect) # plot the points
  }


  if(length(tmp)>1){# if more than one sample= compute the ""spider" around the centroid
    # get coordinates
    co <- coordinates(tmp)
    Me <- colMeans(co) # get the mean locations for the ACS
    AC <- SpatialPoints(cbind(Me[1], Me[2]))

    ## add the segments#
    if(plot.segments==TRUE){
      for(j in 1:length(tmp)){
       segments(x0=Me[1], x1=co[j,1], y0=Me[2], y1=co[j,2], col=col.segments[i] ,lwd=lwd.segm)
      }
    }

    points(tmp, pch=pch.detect, col=col.detections[i],cex=cex, bg=bg.detect) # plot the detections

    if(plot.ACS==TRUE){
    points(AC, pch=24, col="black", bg=col.ACS[i],cex=cex) # plot the estimated ACs
    }
  }
 }#i
}



