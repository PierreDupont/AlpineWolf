#' @title MakeHabitat: a function to create a list of objects to define the habitat
#'
#'
#'
#' @description
#' \code{MakeHabitat} returns a \code{List} of objects necessary to define the habitat in SCR (can be slow for large habitat...)
#'
#' @param poly A \code{SpatialPolygon} object with the study area
#'	@param resolution A \code{Numerical} with the grid cell size of the habitat raster.
#' @param buffer A \code{Numeric} with the size of the buffer area around the focal study area (in meters in poly is UTM).
#' @param polygon.clip A \code{SpatialPolygon} defining suitable habitat.
#' @param plot.check A \code{Logical} to display checking plots (if TRUE).
#' @param fasterize A \code{Logical} if TRUE, the fasterize function from the fasterize library is used. Is not used if CoverToKeepHabitat is used
#' @param CoverToKeepHabitat A \code{Numeric} if.null, all cells intersecting with habitat are considered. if CoverToKeepHabitat is \code{Numeric}, this value (expressed in percentage) is used to retain habitat cells with cover >= CoverToKeepHabitat
#'
#' 
#' @return A \code{List} object with the coordinates and attributes of the habitat
#'  @return habitat.sp: A \code{SpatialPointsDataFrame} object with the coordinates of the center of cells of the habitat
#'  @return habitat.clip.sp: A \code{SpatialPointsDataFrame} object with the coordinates of the center of cells of the habitat cliped to the suitable habitat.
#'  @return habitat.xy: A \code{DataFrame} object with the x and y coordinates of the habitat cells.
#'  @return IDCells.mx: A \code{matrix} with the ID of each cell of the raster
#'  @return habitat.r: A \code{raster} with the suitable habitat (1) and non suitable (0)
#'  @return habitat.mx: A \code{matrix} with the suitable habitat (1) and non suitable (0)
#'  @return resolution: A \code{numeric} with the resolution used for the raster 
#'  @return buffered.habitat.poly: A \code{SpatialPolygons} that includes the buffer. 
#'  @return habitat.poly: A \code{SpatialPolygons} with the original polygon used. 
#'  @return habitat.index: A \code{vector} with ID cell of the raster that fall within the suitable habitat. 



MakeHabitat <- function(poly, 			               ## Polygon of the study area (sp object)
                        resolution, 	               ## Resolution of the grid cell (in meters)
                        buffer = NULL, 					## Add a buffer size (in meters)
                        polygon.clip = NULL, 		   ## Add a polygon that would clip over the buffer (for the sea for example )
                        plot.check = TRUE, 			   ## if TRUE some graphics will show up to check if it looks right )
                        fasterize = FALSE,            ## Increase spee, requires library(fasterize)    
                        CoverToKeepHabitat = NULL)
                        { 
   ## ----- Store the original study area polygon -----
   polygon.orig <- poly
   
   ## ----- Create the buffer zone around the study area -----
   if(!is.null(buffer)){
      buff <- rgeos::gBuffer(poly, width = buffer)
      poly <- raster::union(poly, buff)
      }
   
   ## ----- Remove unsuitable habitat from the study area + buffer-----
   if(!is.null(polygon.clip)){
      poly <- gIntersection(polygon.clip, poly, drop_lower_td = TRUE) 
      }
   
   ## ----- Create Habitat raster (study area + buffer)-----   
   r <- raster(extent(poly))                                         ## Rasterise polygon
   res(r) <- resolution  ## Change resolution
   
   ## ----- Keep habitat based on proportion covered -----
   if(!is.null(CoverToKeepHabitat)){
      
      r2 <- raster(extent(poly))                                         ## Rasterise polygon
      res(r2) <- resolution/10  ## Change resolution
      r <- aggregate(r2, fact=10)
      
      if(fasterize==TRUE){
         polygon.r <- fasterize(st_as_sf(poly), r2)
         polygon.r[is.na(polygon.r)]<-0
         polygon.r <- aggregate( polygon.r, fact=10, sum, na.rm=TRUE)
          }else{
            polygon.r <- rasterize(poly, r, getCover=TRUE)        ## Add field=1 so its 1 for study area and NA for the rest.
          }
      polygon.r[polygon.r < CoverToKeepHabitat] <- 0
      polygon.r[polygon.r >= CoverToKeepHabitat] <- 1
      
   }else{
   ## ----- Use fasterize (a bit faster)-----
    if(fasterize==TRUE){
      polygon.r <- fasterize(st_as_sf(poly), r)
      polygon.r <- rasterize(as(poly, "SpatialLines"), polygon.r,update=T, field=1)      ## rasterise as lines and then polygons so all cells touched by the polygon are covered 
    }else{
       polygon.r <- rasterize(as(poly, "SpatialLines"), r, field=1)      ## rasterise as lines and then polygons so all cells touched by the polygon are covered 
       polygon.r <- rasterize(poly, polygon.r, update=T, field=1)        ## Add field=1 so its 1 for study area and NA for the rest.
    }
   
   }
   
   ## ----- Create Habitat matrix (habitat : 1 and non-habitat: 0) -----
   habitat.r <- polygon.r
   habitat.r[is.na(habitat.r)] <- 0                                     ## Give 0 values to raster cells outside the study area
   habitat.mx <- as.matrix(habitat.r)                                   ## Convert to matrix
   
   ## ----- Give unique IDs to cells ----- 
   IDCells.r <- polygon.r
   IDCells.r[] <- 1:length(IDCells.r)							         ## Cell ID starts from the top left corner 
   IDCells.mx <- as.matrix(IDCells.r)							               ## Convert to matrix
   
   ## ----- Obtain xy coordinates of cells -----   
   habitat.xy <- xyFromCell(polygon.r, 1:ncell(polygon.r))
   dimnames(habitat.xy) <- list(1:length(habitat.xy[,"x"]), c("x","y"))
   habitat.sp <- SpatialPointsDataFrame(data.frame(habitat.xy[,c("x","y")]), data=data.frame(habitat.xy), proj4string=CRS(projection(poly)))
   habitat.index <- which(!is.na(over(as(habitat.sp, "SpatialPoints"), as(poly,"SpatialPolygons"))))
   habitat.clip.sp <- habitat.sp[!is.na(over(as(habitat.sp, "SpatialPoints"), as(poly, "SpatialPolygons"))), ]
   
   ## -- Obtain lower and upper cell coordinates
   lower.hab.sp <- data.frame(coordinates(habitat.sp) - resolution/2)
   upper.hab.sp <- data.frame(coordinates(habitat.sp) + resolution/2)
   coordinates(upper.hab.sp) <- upper.hab.sp
   coordinates(lower.hab.sp) <- lower.hab.sp
   proj4string(lower.hab.sp) <- CRS(projection(poly))
   proj4string(upper.hab.sp) <- CRS(projection(poly))
   
   ## ----- Visual plotting to check if everything is right -----  
   if(plot.check)
      {
      plot(habitat.r)							                          ## Visual check of the habitat and polygons 
      plot(poly, add=TRUE, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))						
      plot(polygon.orig, add=TRUE, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 1))
      }
      
      ## ----- List of output objects -----
   output <- list(habitat.sp = habitat.sp,
                 habitat.clip.sp = habitat.clip.sp,
                 habitat.xy = habitat.xy,
                 IDCells.mx = IDCells.mx,
                 habitat.r = habitat.r,
                 habitat.mx = habitat.mx,
                 resolution = resolution,
                 buffered.habitat.poly = poly,
                 habitat.poly = polygon.orig,
                 habitat.index = habitat.index,
                 upper.hab.sp = upper.hab.sp,
                 lower.hab.sp = lower.hab.sp
                 #----PLEASE DO NOT REMOVE (RB)
                 )
      return(output)
   }
   
   