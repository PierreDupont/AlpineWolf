
#' @title Identify clusterered NGS data in Space and time 
#'
#' @description
#' \code{CleanData} returns a \code{SpatialPointsDataFrame} object of the same dimensions than the input one, but with a column "Delete", to check whether 
#' it should be considered as a duplicate basedeted based on space and time.
#'
#' @param myData.sp A \code{SpatialPointsDataFrame} object containig  DNA samples data, basically after runing the \code{CleanData} function.
#'	@param resolution A \code{Numeric} value for the resolution of the raster .
#' @param time.interval A \code{string} containing the time interval to be used: time interval can also be "sec", "min", "hour", "day", "DSTday", "week", "month", "quarter" or "year". See \code{seq.POSIXt}
#' @param remove.dupli A \code{string} TRUE: if duplicated raw in time and space should be removed from the SpatialPointsDataframe or not (FALSE). 
#' 
#' @return A \code{SpatialPointsDataFrame} object with the x and y locations of the different samples with associated attributes 
#' Id: the individual id
#' DNAID: the sample Id
#' Year: the sample biological year
#' Date: the actual dinfing date  
#' Species: the identified species
#' Country: the country the sample was found in
#' Birth. alleged birth year (when available; based on dead recoveries)
#' Death: biological year of death (for dead recoveries)
#' Age: alleged age (when available; based on dead recoveries)
#' Delete: whether the row is conisdered as a duplicate or not. 


IdentifyClusteredNGS <- function( myData.sp
                                , resolution = 100      # in meters 
                                , time.interval= "day"  # time interval can also be "sec", "min", "hour", "day", "DSTday", "week", "month", "quarter" or "year"
                                , remove.dupli = TRUE)
{
   # get list of IDs
   IDS <- as.character(unique(myData.sp$Id))
   
   # create a raster
   r <- raster((extent(myData.sp)+resolution), resolution=resolution)
   proj4string(r) <- CRS(proj4string(myData.sp))
   values(r) <- 1:ncell(r)
   
   # extract the cellsID for each sample
   myData.sp$cellIDS <- cellFromXY(r, myData.sp) 
   xycell <- xyFromCell(r, myData.sp$cellIDS)
   
   ## add a delete column
   myData.sp$Delete <- "NO"
   myData.sp$X <- xycell[ ,1]
   myData.sp$Y <- xycell[ ,2] 
   # convert as a data.frame
   myData.sp.df <- data.frame(myData.sp)[ ,1:ncol(myData.sp)]
   
   
   # obtain sequence of dates 
   date.range <- range(myData.sp$Date, na.rm = T)
   date.range[2] <-  date.range[2]+ 86400*15 ## add some days so the last date is also covered by the interval 
   seq.date <- seq(date.range[1], date.range[2], by= time.interval)
   
   
   # start for each time interval 
   for(i in 2:length(seq.date)){
      tmp  <- myData.sp[myData.sp$Date >= seq.date[ (i-1)] & myData.sp$Date < seq.date[i], ]
      
      if(length(tmp) > 1 ){
         
         IDS <- unique(tmp$Id)
         
         #for each ids 
         for ( j in 1:length(IDS)){
            
            tmp2 <- tmp[tmp$Id==IDS[j],]
            uni.cell <- unique(tmp2$cellIDS)
            
            if(length(uni.cell) < length(tmp2) ){
               # for each cells #
               tab <- Filter(function(x) x > 1, table(tmp2$cellIDS))
               
               uni.cell1 <- as.numeric(names(tab)[[1]])
               for (k in 1:length(uni.cell1)){
                  tmp3 <- tmp2[tmp2$cellIDS==uni.cell1[k], ]
                  
                  if(length(tmp3) >1 ){
                     myData.sp.df[dimnames(tmp3@data)[[1]][2:length(tmp3)], ]$Delete <- "YES"# add yes if you want remove the duplicates
                  }
               }#k
            }   
         }#j  
      }
      #print(i)
   }#i
   
   # make it sp 
   myData.sp.df.sp <- SpatialPointsDataFrame(myData.sp.df[,c("X","Y")], myData.sp.df[ ,(1: (ncol(myData.sp.df)-2)) ]  )
   proj4string(myData.sp.df.sp) <- CRS(proj4string(myData.sp))
   # if duplicates should be removed 
   if(remove.dupli==TRUE){
      print(paste(nrow(myData.sp.df.sp[myData.sp.df.sp$Delete=="YES",]), "detections have been removed",sep=" "))
      myData.sp.df.sp <-   myData.sp.df.sp[myData.sp.df.sp$Delete=="NO",]
      
   }
   
   return(myData.sp.df.sp)
}
