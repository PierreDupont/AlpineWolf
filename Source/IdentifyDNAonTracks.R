#' @title Function to identify DNA samples found on SNO search tracks 
#' @description
#' \code{IdentifyDNAonTracks} returns a list with: 
#'                          I.samples found on search tracks,
#'                          II.sample found outside search tracks,
#'                          III. the percentage of DNA samples found on these tracks
#' 
#' @param myData A \code{SpatialLinesDataFrame} sp class object with DNA samples' locations
#' @param myTracks A \code{SpatialLinesDataFrame} sp class object with GPS search tracks
#' @param resolution A \code{numeric} to set the maximum distance between a sample and a search track 
#' @param exact.date A \code{logical} determines if DNA samples should match the exact dates of search tracks 
#' (default is FALSE meaning DNA samples may be assigned to tracks even if they were found before or after the actual date of the track)
#' @param plot.check A \code{logical} to display validation plots.
#' @param map	A \code{polygon} to plot the study area considered.
#' @keywords SCR data prep
#'
#' @return A \code{list} object containing points on tracks, points outside tracks and proportion of points on the tracks.

IdentifyDNAonTracks <- function( myData 
                               , myTracks 
                               , resolution = 100
                               , exact.date = FALSE
                               , plot.check = TRUE
                               , map = map)
  {
   ##-----------------------------------------------------------------------------------------------
   if(!exact.date){
      myDist <- NULL
      for(i in 1:dim(myData)[1]){
         myDist <- c(myDist,gDistance(myData[i,], myTracks))
         }#i                                   
      ScatsOnTracks <- myData[myDist <= resolution, ] 
      if(plot.check){
         myCol <- ifelse(myDist <= resolution, "green","red")
         plot(map)                                             # shows study area 
         plot(myTracks, col="gray60", add=TRUE)                # shows search tracks 
         plot(myData, pch = 21, bg = myCol, add=TRUE)          # shows all DNA samples
         plot(map, add = TRUE)
         }#if
      }#if
  
   ##-----------------------------------------------------------------------------------------------
   ### If we want the exact dates to match not just the years
   if(exact.date){
      myTracks$Date <- as.POSIXct(strptime(myTracks$DATE_DEF, "%Y/%m/%d"))
      temp.data <- myData[myData$Date %in% myTracks$Date, ]
      myDist <- myCol <- NULL
      if(plot.check){ 
         plot(map)                                                                        
         plot(myTracks, col="gray60", add=TRUE)
         }#if
      if(dim(temp.data)[1]>=1){
         for(i in 1:dim(temp.data)[1]){
            temp.tracks <- myTracks[myTracks$Date == temp.data$Date[i], ]
            myDist <- c(myDist, gDistance(temp.data[i,],temp.tracks))
            myCol <- c(myCol, ifelse(myDist[i] <= resolution,"green","red"))
            }#i
         }#if
      ScatsOnTracks <- temp.data[myDist <= resolution, ]   
      plot(temp.data, add=TRUE, pch=21, bg = myCol)
      }
   
   ##--------------------------------------------------------------------------------------------------
   ## OutPut files 
   ScatsOutTracks <- myData[!myData$DNAID %in% ScatsOnTracks$DNAID, ]# DNA samples NOT on search tracks
   Prop <- dim(ScatsOnTracks)[1]/dim(myData)[1]                      # Proportion of points ON tracks
   
   return(list( ScatsOnTracks.sp = ScatsOnTracks
              , ScatsOutTracks.sp = ScatsOutTracks
              , Proportion = Prop))
   }