#' @title Generate Detection Histories Arrays
#'
#' @description
#' \code{MakeY} returns a \code{Array} object with the detections (presence/absence or counts) of individuals over detectors and time.
#'
#' @param myData A \code{SpatialPointsDataFrame} of detections with ID, Year and Detector columns.
#'	@param myDetectors A \code{list} of SpatialPointsDataFrames containing the detector locations.
#' @param method A \code{Character} object defining the array to be returned ("Binomial" , "Poisson" , "Bernoulli").
#' @param myData2 OPTIONAL: A second \code{SpatialPointsDataFrame} of detections with ID, Year and Detector columns.
#'	@param myDetectors2 OPTIONAL: A \code{SpatialPointsDataFrame}of SpatialPointsDataFrames containing the detector locations.
#' @param method2 OPTIONAL: A \code{Character} object defining the array to be returned ("Binomial" , "Poisson" , "Bernoulli").
#' @param Id OPTIONAL: A \code{vector} of IDs to return the y array for.
#' @param returnIdvector OPTIONAL: A \code{logical} to return the vector of individual ID.
#' 
#' @return A \code{Array} object with the individual in rows, the detectors in columns and years in third dimension. 

MakeY <- function( myData,                 
                   myDetectors,           
                   method = "Bernoulli",
                   myData2 = NULL, 
                   myDetectors2 = NULL,
                   method2 = "Bernoulli",
                   Id = NULL,
                   returnIdvector = FALSE){
   
   ## PREPARE THE DATA 
   myData <- as.data.frame(myData)
   myData2 <- as.data.frame(myData2)
   
   if(sum(duplicated(myData2$Id))>0){stop("INDIVIDUALS CANNOT BE DEAD TWICE!!!!")}
   
   if(!is.list(myDetectors)){myDetectors <- list(myDetectors)}
   if(!is.list(myDetectors2)){myDetectors2 <- list(myDetectors2)}
   
   ## AGGREGATE MULTIPLE DETECTIONS AT THE SAME SUBDETECTOR TO ONE DETECTIONS ONLY
   if(method == "Binomial"){
      myData <- myData[!duplicated(myData[ ,c("sub.detector","Year", "Id" )]), ] 
      myData <- droplevels(myData)
   }
   
   ## NUMBER OF INDIVIDUALS
   if(is.null(Id)){
      if(length(myData2)>0){
         Id <- unique(c(as.character(myData$Id), as.character(myData2$Id)))
      } else {Id <- unique(as.character(myData$Id))}
   }
   N_Id <- length(Id)
   
   ## NUMBER OF YEARS
   Years <- sort(unique(c(myData$Year, myData2$Year)))#[RB] ADDED SORTING 2019-07-30 - reason: ran into problems with bear data
   N_Years <- length(Years)
   
   ## NUMBER OF DETECTORS 
   N_Detec <- max(unlist(lapply(myDetectors, function(x) dim(x)[1])))
   Detectors <- 1:N_Detec
   
   
   ##--[RB] FASTER/LESS HUNGRY DUMMY TABLE CREATION
   
   
   dummy<-myData[ ,c("Id","Year","Detector")]
   dummy$Id<-factor(dummy$Id,Id)
   dummy$Year<-factor(dummy$Year,Years)
   dummy$Detector<-factor(dummy$Detector,Detectors)
   temp <- table(dummy$Id, dummy$Detector, dummy$Year)
   #dim(temp)
   
   # ## CREATE A DUMMY DATASET
   # D_Id <- rep(Id, N_Years*N_Detec)
   # D_Year <- rep(rep.int(Years, rep(N_Id,N_Years)), N_Detec)
   # D_Detector <- rep.int(Detectors, rep(N_Id*N_Years, N_Detec))
   # dummy <- cbind.data.frame(Id = D_Id, Year = D_Year, Detector = D_Detector)
   # 
   # ## COMBINE DUMMY AND REAL DATA
   # dummy <- rbind(dummy, myData[ ,c("Id","Year","Detector")])
   # 
   # ## CREATE THE DETECTION ARRAY
   # temp <- table(dummy$Id, dummy$Detector, dummy$Year)
   id_names <- row.names(temp)
   Count.ar <- array(NA, c(dim(temp)))
   Count.ar[] <- temp#-1
   dimnames(Count.ar) <- dimnames(temp)
   if(method =="Bernoulli"){Count.ar[Count.ar>=1] <- 1}
   
   
   # [CM] REORDER BY INDIVIDUALS TO MATCH THE PREVIOUS VERSION OF MAKE Y 
   Count.ar <- Count.ar[order(dimnames(Count.ar)[[1]]),,]
   id_names <- id_names[order(id_names)]
   
   ## OUTPUT SETTINGS
   if(returnIdvector){
      output <- list(y.ar = Count.ar, Id.vector = id_names)
   } else {output <- list(y.ar = Count.ar)}
   
   ## IF SECOND DATASET (e.g. DEAD RECOVERIES)
   if(length(myData2)>0){
      ## AGGREGATE MULTIPLE DETECTIONS AT THE SAME SUBDETECTOR TO ONE DETECTIONS ONLY
      if(method2 == "Binomial"){
         myData2 <- myData2[!duplicated(myData2[ ,c("sub.detector","Year", "Id" )]), ] 
         myData2 <- droplevels(myData2)
      }
      
      # NUMBER OF DETECTORS 
      N_Detec2 <- max(unlist(lapply(myDetectors2, function(x) length(x))))
      Detectors2 <- 1:N_Detec2
      
      ##--[RB] FASTER/LESS MEMORY-HUNGRY DUMMY TABLE CREATION

      dummy2<-myData2[ ,c("Id","Year","Detector")]
      dummy2$Id<-factor(dummy2$Id,Id)
      dummy2$Year<-factor(dummy2$Year,Years)
      dummy2$Detector<-factor(dummy2$Detector,Detectors)
      temp2 <- table(dummy2$Id, dummy2$Detector, dummy2$Year)
      #dim(temp2)
      
      
      # ## CREATE A DUMMY DATASET
      # D_Id2 <- rep(Id, N_Years*N_Detec2)
      # D_Detector2 <- rep.int(Detectors2, rep(N_Id*N_Years, N_Detec2))
      # D_Year2 <- rep(rep.int(Years, rep(N_Id,N_Years)), N_Detec2)
      # dummy2 <- cbind.data.frame(Id = D_Id2, Year = D_Year2, Detector = D_Detector2)
      # 
      # ## COMBINE DUMMY AND REAL DATA
      # dummy2 <- rbind(dummy2, myData2[ ,c("Id","Year","Detector")])
      # 
      # ## CREATE THE DETECTION ARRAY
      # temp2 <- table(dummy2$Id, dummy2$Detector, dummy2$Year)
      Count.ar2 <- array(NA,c(dim(temp2)))
      Count.ar2[] <- temp2#-1
      dimnames(Count.ar2) <- dimnames(temp2)
      if(method=="Bernoulli"){Count.ar2[Count.ar2 >= 1] <- 1}
      
      
      
      # ## OUTPUT SETTINGS
      # if(returnIdvector){
      #    output <- list(y.ar = as.array(Count.ar), y.ar2 = as.array(Count.ar2), Id.vector = id_names)
      # } else {output <- list(y.ar = Count.ar, y.ar2 = Count.ar2)}
      # 
      
      ## OUTPUT SETTINGS
      Count.ar <- as.array(Count.ar)
      Count.ar2 <- as.array(Count.ar2)
      # reorder to match previous MAKE Y 
      Count.ar <- Count.ar[order(dimnames(Count.ar)[[1]]),,]
      Count.ar2 <- Count.ar2[order(dimnames(Count.ar2)[[1]]),,]
      id_names <- id_names[order(id_names)]
      
      if(returnIdvector){
         output <- list(y.ar = as.array(Count.ar), y.ar2 = as.array(Count.ar2), Id.vector = id_names)
      } else {output <- list(y.ar = Count.ar, y.ar2 = Count.ar2)}
      
      
      
      
   }#if
   return(output)
}

