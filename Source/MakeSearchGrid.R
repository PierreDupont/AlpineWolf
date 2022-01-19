#' @title Function to create detectors
#' #'
#' @description
#' \code{MakeSearchGrid} creates  detectors within a defined area (data) and resolution between detectors. 
#' Division precises if subspatial division should be performed (Following PAB model).  
#' 
#' @param data A \code{sp} object. can be a \code{SpatialPolygons} or \code{SpatialPoints} or \code{SpatialLines} or \code{raster}. Spatial Lines objects takes time to compute
#' @param resolution Numeric variable denoting the size of grid cells in units of (\code{polygon}).
#' @param div Numeric variable denoting the number of equally-sized subplots grid cell is to be divided into.
#' @param center Logical variable denoting whether nodes of the resulting search grid are to be centered within their respective subplots.
#' @param plot Logical for whether \code{True} or not (\code{FALSE}) check plots are to be generated.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/MakeSearchGrid.R
#' @keywords simul
#'
#' @examples
#' # Make a nested search grid:
#' MakeSearchGrid(...)

MakeSearchGrid <- function(data = data,
                           resolution = resolution,
                           div = 1,
                           center = T,
                           plot = TRUE,
                           fasterize = FALSE){
   
   ## if the data is already a raster   
   subdetector.r <- data
   
   ## Obtain the resolution of the subdetectors 
   res1 <- resolution/sqrt(div)  
   
   ## Convert other data types to easier simpler format (no dataframe). 
   if(class(data)=="SpatialPolygonsDataFrame"){data <- as(data,"SpatialPolygons")}
   if(class(data)=="SpatialPointsDataFrame"){data <- as(data,"SpatialPoints")}
   if(class(data)=="SpatialLinesDataFrame"){data <- as(data,"SpatialLines")}
   
   ### ==== CREATE THE SUBDETECTORS  ====
   ## if SPATIALLINES
   if(class(data)=="SpatialLines"){
      ## Alternative of using rasterize. should be a bit faster for lines 
      detector.r1 <- raster(extent(data), resolution=res1, crs=proj4string(data) )
      detector.r1[] <- 1:ncell(detector.r1)
      
      ## Intersect lines with raster "polygons" and add length to new lines segments
      rsp <- rasterToPolygons(detector.r1)
      rp <- raster::intersect(data, rsp)
      rp$length <- gLength(rp, byid=TRUE)
      x <- tapply(rp$length, rp$layer, sum)
      subdetector.r <- raster(detector.r1)
      subdetector.r[as.integer(names(x))] <- x
   }
   
   ## if SPATIALPOINTS 
   if(class(data)=="SpatialPoints"){
      detector.r1 <- raster(extent(data), resolution = res1, crs = proj4string(data))
      subdetector.r <- rasterize(data, detector.r1, fun = "count")
   }
   
   ## if SPATIALPOLYGONS
   if(class(data)=="SpatialPolygons"){
      detector.r1 <- raster(extent(data), resolution = res1, crs = proj4string(data))
      
      if(fasterize==TRUE){
         subdetector.r <- fasterize(st_as_sf(data), detector.r1)
      }else{
         subdetector.r <- rasterize(data, detector.r1)
      }
   }
   
   ### ==== AGGREGATE TO MAIN DETECTORS  ====
   if(sqrt(div)>1){
      maindetector.r <- aggregate(subdetector.r, fact = sqrt(div), fun = sum)
   }else{
      maindetector.r <- subdetector.r
   }
   
   ### ==== CREATE POLYGONS FROM RASTER  ====
   ## Main detectors 
   temp.r <- raster(maindetector.r)
   temp.r[] <- 1:length(temp.r)
   maindetector.poly <- rasterToPolygons(temp.r, dissolve = TRUE)
   
   ## Sub-detectors
   temp.r <- subdetector.r
   temp.r[] <- 1:length(temp.r)
   subdetector.poly <- rasterToPolygons(temp.r, dissolve = TRUE)
   
   ### ==== OBTAIN SPATIALPOINTS FROM DETECTORS ====
   ## Main detectors 
   main.detector.xy <- xyFromCell(maindetector.r, 1:ncell(maindetector.r))
   main.detector.sp <- SpatialPointsDataFrame( data.frame(main.detector.xy[,c("x","y")]),
                                               data=data.frame(main.detector.xy),
                                               proj4string=CRS(projection(data)))
   names(main.detector.sp@data) <- c("main.cell.x","main.cell.y")
   main.detector.sp@data$main.cell.id <- 1:length(main.detector.sp)
   
   ## Sub-detectors 
   if(center){
      ## ALTERNATIVE 1: center points of subdetector detector cells
      detector.xy <- xyFromCell(subdetector.r, 1:ncell(subdetector.r))
      sub.detector.sp <- SpatialPointsDataFrame( data.frame(detector.xy[ ,c("x","y")]),
                                                 data = data.frame(detector.xy),
                                                 proj4string = CRS(projection(data)))
   }else{
      ## ALTERNATIVE 2: one random point within each subdector cell
      sub.detector.sp <- lapply(1:length(subdetector.poly), function(x){
         spsample(subdetector.poly[x, ], 1, type="random")
      })
      sub.detector.sp <- do.call(rbind, sub.detector.sp)  
      sub.detector.sp <- SpatialPointsDataFrame( coordinates(sub.detector.sp),
                                                 data = data.frame(coordinates(sub.detector.sp)),
                                                 proj4string = CRS(projection(data)))
   }
   
   names(sub.detector.sp@data) <- c("x","y")
   sub.detector.sp@data$Id <- 1:length(sub.detector.sp)
   
   ### ==== SELECT ACTIVE DECTECTORS (SEARCHED) ====
   ## merge information so we get id and x and y coordinates of main detectors in the subdetector file
   sub.detector.sp$main.cell.id <- over(sub.detector.sp, maindetector.poly)[,1]
   merge.df <- merge(sub.detector.sp@data, main.detector.sp@data, by="main.cell.id")
   merge.df <- merge.df[order(merge.df$Id), ]
   sub.detector.sp@data <- merge.df
   
   ## Subset main detectors 
   values.main.r <- raster::extract(maindetector.r, main.detector.sp)
   values.main.r[is.na(values.main.r)] <- 0
   main.detector.sp <- main.detector.sp[values.main.r>0, ]
   main.detector.sp$count <- values.main.r[values.main.r>0]
   
   ## Subset sub detectors 
   values.sub.r <- raster::extract(subdetector.r, sub.detector.sp)
   values.sub.r[is.na(values.sub.r)] <- 0
   sub.detector.sp <- sub.detector.sp[values.sub.r>0, ]
   sub.detector.sp$count <- values.sub.r[values.sub.r>0]
   
   ### ==== MAKE A PLOTTING FUNCTION ====
   if(plot){
      if(sum(class(data) %in% c("SpatialPolygons","SpatialPoints","SpatialLines"))>0){
         plot(data, col=grey(0.3))
         ## else: raster layer
      }else{plot(data)}
      plot(maindetector.poly, add=TRUE, lwd=3)
      plot(subdetector.poly, add=TRUE)
      plot(sub.detector.sp, pch=19, cex=0.6, col=as.numeric(sub.detector.sp@data$main.cell.id), add=TRUE)
      points(main.cell.y ~ main.cell.x, data=sub.detector.sp@data, pch=19, cex=1, col=as.numeric(sub.detector.sp@data$main.cell.id))
   }
   
   ### ==== OUTPUT ====
   out <- list( detector.sp = sub.detector.sp,
                main.detector.sp = main.detector.sp,
                sub.grid.poly = subdetector.poly,
                grid.poly = maindetector.poly)
}


