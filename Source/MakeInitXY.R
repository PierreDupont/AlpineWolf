#' @title Function to set initial values of XY coordinates of ACS
#'
#' @description
#' \code{InitXY} returns a matrix object with with the x coordinates ([,1]) and y coordinates  ([,2]). it returns the location of a detection for an detected individual and a random location for augmented individuals. 
#' 
#' @param y \code{matrix} or \code{array}  with individual detections. row=Ids, col= detectors and if arrray=  [,,t] time. 
#' @param habitat.mx \code{matrix} from MakeHabitat
#' @param detector.xy \code{matrix} with coordinates of detectors scaled
#' @param IDCells.mx \code{matrix} the ID of the cells of the habitat
#' @param grid.xy \code{matrix} with the scaled coordinates of the habitat 
#' @param xy.bounds \code{matrix} with the scaled xy bounds of the moving window when Fragmented SCR is used.
#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' MakeInitXY(   y = data$y
#'              , habitat.mx = data$habitat.mx
#'              , detector.xy = data$detector.xy
#'              , grid.xy = grid.xy
#'              , xy.bounds=NULL )

MakeInitXY <- function(    y = y
                         , habitat.mx = habitat.mx
                         , detector.xy = detector.xy
                         , IDCells.mx = IDCells.mx
                         , grid.xy = grid.xy
                         , xy.bounds = NULL ){
   
   # ---- STEP 1: GET THE OBJECT READY TO STORE THE DATA ----- 
   # Make it general to work with the time dimention 
   n.years <- ifelse(length(dim(y))>=3, dim(y)[3], 1)
   n.individuals <- dim(y)[1]
   n.detectors <- dim(detector.xy)[1]

   if(length(dim(y)) == 2){y <- array(y, c(n.individuals, n.detectors, n.years))}
   if(length(dim(detector.xy)) == 2){detector.xy <- array(detector.xy, c(n.detectors, 2, n.years))}
   n.year <- dim(y)[3]
   
   
   # STORE THE DIMENSION OF THE WINDOW 
   # IF NO FRAGMENTATION IS USED 
   dim.mx <- array(0, c(n.individuals, 2, 2, n.year))
      for(t in 1:n.year){
      dim.mx[,1,1,t] <- 1 # min x
      dim.mx[,2,1,t] <- 1 # min y
      dim.mx[,1,2,t] <- dim(habitat.mx)[2]# max x
      dim.mx[,2,2,t] <- dim(habitat.mx)[1]# max y 
         
         if(!is.null(xy.bounds)){
             dim.mx <- xy.bounds
             dim.mx[,1,1,t] <-  floor(dim.mx[,1,1,t]) + 1 # Add one to match the habitat 
             dim.mx[,2,1,t] <-  floor(dim.mx[,2,1,t]) + 1 # Add one to match the habitat
         }
      }
   
   # IF  FRAGMENTATION IS USED 

   
  
   empty <- array(numeric(),c(n.individuals, 2, n.year))
   ids <- lapply(c(1:n.individuals), function(x) x)
   
   # ---- STEP 2: FIND SXY FOR ALL INDIVIDUALS  ----- 
   
      for(t in 1:n.year){
         listt <- apply(y[,,t], 1, function(x) which(x>0, arr.ind = T))# obtain detectors with detections for each ids 
        
          if(length(listt)==0){listt <- list()
            for(i in 1:n.individuals){
                listt[[i]] <- integer()
             }
          }
         
          empty[,,t] <- do.call(rbind, lapply(ids, function(i) {
             
            #print(i)
            x <-  listt[[i]]
            # if only one detections 
            if(length(x)==1){ # If only one detection, use that detection as a starting value
               detector.xy[x,,t]
            }else{
            
            # if several detections    
            if(length(x)>1){# if more than 1 detection, use average value
                  mean.det <- colMeans(detector.xy[x,,t])
                  # if doesnt end up in habitat, use a random coordinate within habitat 
                  if(habitat.mx[floor(mean.det)[2]+1,floor(mean.det)[1]+1]==1){
                    mean.det 
                  }else{
                     min.x <- dim.mx[i,1,1,1]
                     max.x <- dim.mx[i,1,2,1]
                     min.y <- dim.mx[i,2,1,1]
                     max.y <- dim.mx[i,2,2,1]
                     
                     window.mx.ind <- habitat.mx[c(min.y:max.y), c(min.x:max.x)]
                     window.ID.ind <- IDCells.mx[c(min.y:max.y), c(min.x:max.x)]
                     
                     id.cell <- window.ID.ind[window.mx.ind==1]
                     sxy.grid <- matrix(grid.xy[id.cell,], nrow=length(id.cell),ncol=2)
                     
                     matrix(sxy.grid[floor(runif(1, 1, length(id.cell))),], ncol=2, nrow=1)                     }
            }else{
            
            # if no detections (augmented IDS)
            if(length(x)==0){
               
               min.x <- dim.mx[i,1,1,1]
               max.x <- dim.mx[i,1,2,1]
               min.y <- dim.mx[i,2,1,1]
               max.y <- dim.mx[i,2,2,1]
               
               window.mx.ind <- habitat.mx[c(min.y:max.y), c(min.x:max.x)]
               window.ID.ind <- IDCells.mx[c(min.y:max.y), c(min.x:max.x)]
                                
               id.cell <- window.ID.ind[window.mx.ind==1]
               sxy.grid <- matrix(grid.xy[id.cell,], nrow=length(id.cell),ncol=2)
               
               matrix(sxy.grid[floor(runif(1, 1, length(id.cell))),], ncol=2, nrow=1)
                }
             }}      
          })) 
      
      }


   
   
   return(empty)
   
}



