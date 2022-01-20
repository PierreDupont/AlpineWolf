#' @title Function to extract the most likely AC for each individual 
#' #'
#' @description
#' \code{GetMostLikelyAC} returns a dataframe and sp object with Mostlikely AC for each ids
#' 
#' @param ID.m \code{matrix} with the cell IDs created from MakeHabitat
#' @param posterior.sy \code{SpatialPointsDataframe} Posterior sy of  jags. is typically a matrix with individuals (i) in columns i.e. sy[,i] 
#' @param posterior.sx \code{SpatialPointsDataframe} Posterior sx of  jags. is typically a matrix with individuals (i) in columns i.e. sx[,i] 
#' @param posterior.z \code{SpatialPointsDataframe} Posterior z of  jags. is typically a matrix with individuals (i) in columns i.e. z[,i] 
#' @param study.area.xy \code{matrix} with the coordinates of the habitat in its original coordinates
#' @param study.area.scaled.xy \code{matrix} with the coordinates of the habitat in its scaled coordinates
#'
#' @examples
#' # UTMToGrid:
#' 
#'GetMostLikelyAC(ID.m = ID.m, posterior.sy = posterior.sy, posterior.sx = posterior.sx, posterior.z = posterior.z,
#'                study.area.xy = study.area.xy, study.area.scaled.xy = study.area.scaled.xy,IDs= c(1:10) )
#
#' 
#' 



GetMostLikelyAC <- function(
  ID.m = ID.m,
  posterior.sy = posterior.sy, #  posteriors sy of  jags. is typically a matrix with individuals (i) in columns i.e. sy[,i] 
  posterior.sx = posterior.sx, #  posteriors sx of jags. is typically a matrix with individuals (i) in columns i.e. sx[,i] 
  posterior.z = posterior.z, #  posteriors z of jags. is typically a matrix with individuals (i) in columns i.e. z[,i] 
  study.area.xy = study.area.xy,
  study.area.scaled.xy = study.area.scaled.xy,
  IDs= c(1:10) 
){
  
  df <- data.frame(ID = 0, cellID=0, x.scaled=0, y.scaled=0, x.origin=0, y.origin=0 )
  df.list <- list()
  
  for(i in 1:length(IDs)){
     
     cell_ID <- 0
     
    
    for ( j in 1:length(posterior.sy[ ,IDs[i]])){
       if(posterior.z[j,IDs[i]] ==1){
      coo <- c(trunc(posterior.sy[j, IDs[i]]+1), trunc(posterior.sx[j, IDs[i]]+1))
      cell_ID[j] <- ID.m[coo[1], coo[2]]
       }
    }#j
    
    df$ID <- IDs[i]
    df$cellID <- as.numeric(dimnames(sort(table(cell_ID), decreasing = T))$cell_ID[1])
    df$x.scaled <- study.area.scaled.xy[df$cellID, 2]
    df$y.scaled <- study.area.scaled.xy[df$cellID, 1]
    df$x.origin <- study.area.xy[df$cellID, 2]
    df$y.origin <- study.area.xy[df$cellID, 1]
    df$z <- mean(posterior.z[, IDs[i]])
    df.list[[i]] <- df
  }#i
  
  
  
  results <- do.call(rbind, df.list)
  
  list.return <- list(predictedAC.xy = results,
                      predictedAC.sp = SpatialPointsDataFrame(results[,c("x.origin","y.origin")], results) )
    
  
  
  return(list.return)
}
