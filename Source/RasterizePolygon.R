RasterizePolygon<-function(poly,r,CoverToKeepHabitat=0.5,fasterize=FALSE){
   
   ## ----- Create Habitat raster (study area + buffer)-----   
   #r <- raster(extent(poly))                                         ## Rasterise polygon
   #res(r) <- resolution  ## Change resolution
   
   r<-raster::trim(r)
   ## ----- Keep habitat based on proportion covered -----
   if(!is.null(CoverToKeepHabitat)){
      
      r [r!=1]<-NA
      r2 <- r#raster(extent(poly))                                         ## Rasterise polygon
      res(r2) <- res(r2)/10  ## Change resolution


      if(fasterize==TRUE){
         polygon.r <- fasterize(st_as_sf(poly), r2)
         polygon.r[is.na(polygon.r)]<-0
         polygon.r <- aggregate( polygon.r, fact=10, sum, na.rm=TRUE)
         #new.r<-overlay(r,polygon.r,fun=sum)
         polygon.r[]<-ifelse(!is.na(r[]),polygon.r[],NA)
      }else{
         polygon.r <- rasterize(poly, r, getCover=TRUE)        ## Add field=1 so its 1 for study area and NA for the rest.
      }
      

      polygon.r[polygon.r < CoverToKeepHabitat] <- NA
      polygon.r[polygon.r >= CoverToKeepHabitat] <- 1
      #plot(polygon.r)
      
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
   
   return(polygon.r)
   
}