#' @title Function to make a toggle detector grid from a rasterstack and a list of conditions 
#' #'
#' @description
#' \code{MakeToggleDetectors} returns a \code{raster}  object with 1 (active, i.e. searched) and 0 (inactive) cells
#' 
#' @param rasterstack A \code{rasterstack} object with covariates
#' @param resolution A \code{Numeric} denoting whether raster stack should be resample to certain resolution 
#'                  (ngb method is used, see ?resample from raster). if NUll no resampling is performed and original resolution is retained.
#' @param conditions A \code{list} of conditions to turn the grid as active. it should contain the name of the layer as names 
#'                   and the conditions as string i.e. (  list(layer.1 = ">0", layer.2= ">50" ) ) 
#' @param plot.check Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated during simulations.
#' 
#'
#' @examples
#' 
#'  my_toggle <- MakeToggleDetectors(  rasterstack = rasterstack
#'                                   , conditions =  list(  layer.1 = ">0"
#'                                                          , layer.2 = ">50"  )
#'                                   , resolution = NULL
#'                                   , plot.check = TRUE
#'                                   
#'                                   )
#'
#'plot(my_toggle)
#' 
#' 


MakeToggleDetectors <- function(  rasterstack = rasterstack
                                , conditions = conditions
                                , resolution = NULL
                                , plot.check = TRUE
    
 ){

   # create a new empty rasterstack
rstack1 <- rasterstack
rstack1[] <- 0

 # Eval the conditions 
for(i in 1:length(conditions)){
   eval(parse(text=paste("rstack1[[i]][rasterstack[[names(conditions)[i]]]", conditions[[i]],"] <- 1",sep="" ))) 
}

# sum the active grids and put it to 1

sumrstack <- sum(rstack1)
sumrstack[sumrstack > 0] <- 1

# if resample 
if(!is.null(resolution)){
   r <- raster(extent(sumrstack), res=resolution, crs= proj4string(sumrstack))
   sumrstack <- raster::resample(sumrstack, r, method="ngb")
}

# if plot# 
if(plot.check==TRUE){
   plot(sumrstack)
}

return(sumrstack)

}


