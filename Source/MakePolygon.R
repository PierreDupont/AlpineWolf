MakePolygon <- function( coord.x = NULL
                       , coord.y = NULL
                       , max.x = NULL
                       , max.y = NULL
                       , random = FALSE
                       , origin.xy = c(0,0)
                       , n.polygons = NULL
                       , intersect = FALSE
                       , CRSStr = NULL)
{
if(is.null(n.polygons)){n.polygons <- 1}

Srs1 <- list()
for(t in 1:n.polygons){
   if(!is.null(coord.x)){Sr1 <- Polygon(cbind(c(coord.x), c(coord.y)))}

   if(!is.null(max.x) & random){
         x1 <- runif(1, 0, max.x/2)
         y1 <- runif(1, 0, max.y/2)
      
         x2 <- runif(1, max.x/2, max.x)
         y2 <- runif(1, 0, max.y/2)
      
         x3 <- runif(1, max.x/2, max.x)
         y3 <- runif(1, max.y/2, max.y)
      
         x4 <- runif(1, 0, max.x/2)
         y4 <- runif(1, max.y/2, max.y)
      
         coord.x <- c(x1,x2,x3,x4)
         coord.y <- c(y1,y2,y3,y4)
         Sr1 <- Polygon(cbind(c(coord.x), c(coord.y)))
         }
   
   if(!is.null(max.x) & !random){Sr1 <- Polygon(cbind( c(origin.xy[1], origin.xy[1], max.x, max.x)
                                                     , c(origin.xy[2], max.y, max.y, origin.xy[2])))}

   Srs1[[t]] <- Polygons(srl = list(Sr1), ID = paste("StudyArea",t))
   }#t

# create spatial polygons object from lists
if(is.null(CRSStr)){CRSStr <- "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"}
SA <- SpatialPolygons(Srs1, proj4string = CRS(CRSStr))
if(intersect){SA <- rgeos::gIntersection(SA, SA)}
return(SA)
}