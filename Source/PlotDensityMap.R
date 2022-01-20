PlotDensityMap <- function( stack
                             , poly = NULL
                             , titles = NULL
                             , interpolate = FALSE
                             , clip = FALSE
                              , return.raster=FALSE,
                            plot.check=TRUE)
{
   min.col <- min(raster::values(stack))
   max.col <- max(raster::values(stack))
   n.years <- nlayers(stack)
   if(class(stack)[1]=="RasterLayer"){ stack <- stack(stack,stack)}## add this when non-multiple years 
   
   if(clip==TRUE){ 
      stack <- mask(stack, poly )}
   
   
   #par(mfrow=c(ceiling(sqrt(n.years)),ceiling(sqrt(n.years))))
   if(plot.check){
   for (t in 1:n.years) 
   {
      plot( stack[[t]]
            , useRaster = TRUE
            , interpolate = interpolate
            , col = rev(terrain.colors(99))
            , breaks = seq(min.col,max.col,length.out=100)
            , legend = FALSE
            , axes = FALSE)
      if(!is.null(poly)){plot(poly,add=TRUE)}
      plot( raster(stack[,,t])
            , legend.only = TRUE
            , col = rev(terrain.colors(99))
            , breaks = seq(min.col,max.col,length.out=100)
            , legend.width = 1
            , legend.shrink = 0.75
            , axis.args = list( at = seq(round(min.col,1), round(max.col,1), length.out=10)
                                , labels = seq(round(min.col,1), round(max.col,1), length.out=10)
                                , cex.axis=0.6)
            , legend.args = list( text = 'Density'
                                  , side = 4
                                  , font = 2
                                  , line = 2.5
                                  , cex = 0.8))
      if(!is.null(titles)){title(main = titles[t], font.main = 4)}
   }
   }
   if(return.raster)return(stack)
}
