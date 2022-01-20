





RescalePosteriorSxy <- function(  sxy = sxy
                                , original.habitat.r = original.habitat.r
                                , new.raster.resolution = new.raster.resolution
                                , polygon = polygon
                                , dim.xy = NULL # assume it is on the 3rd dimension of the sxy
){
require(fasterize)  
# CREATE A POLYGON OF THE OLD HABITAT
polygon.hab <- st_as_sf(aggregate(rasterToPolygons(original.habitat.r, fun = function(x) x>0)))

# CREATE A RASTER WITH THE NEW RESOLUTION 
new.raster <- raster(extent(polygon.hab))
res(new.raster) <- new.raster.resolution
new.raster[] <- 1

# NOW FASTERIZE THE POLYGON INTO THE NEW RASTER 
new.raster <- fasterize(polygon.hab, new.raster)

# 
matrix.input <- as.matrix(new.raster)
IDCells.r <- new.raster
IDCells.r[] <- 1:length(IDCells.r)							         ## Cell ID starts from the top left corner 
IDCells.mx <- as.matrix(IDCells.r)  

sxy <- GridToUTM(  data.sp = NULL
                     , 
                     grid.sp = SpatialPoints(coordinates(original.habitat.r))
                     , 
                     data.sxy= sxy
                     ,
                     plot.check = F)$data.scaled.xy
  
sxy <- UTMToGrid(  data.sp = NULL
                     , grid.sp = SpatialPoints(coordinates(new.raster))
                     , data.sxy= sxy
                     , plot.check = F)$data.scaled.xy
return(list( sxy = sxy
            , new.raster = new.raster
            , IDCells.mx = IDCells.mx
            , habitat.mx = matrix.input
       ))

}
  
  
  
  
  
                        # plot(polygon)
                       
                        
                       