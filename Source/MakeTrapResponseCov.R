MakeTrapResponseCov <- function( data,
                                 data.dead = NULL,
                                 IDs = NULL){

   ## LIST ALL INDIVIDUALS
   if(is.null(IDs)){
      dummyId <- unique(as.character(data$Id))
      if(!is.null(data.dead)){
         dummyId <- unique(c(as.character(data$Id), as.character(data.dead$Id)))
      }
   } else {
      dummyId <- IDs
      data <- data[data$Id %in% dummyId, ]
      data$Id <- factor(as.character(data$Id), levels = unique(as.character(data$Id)))
   }
   N_dummyId <- length(dummyId)
   
   ## LIST ALL YEARS
   Years <- range(unique(data$Year))
   dummyYears <- Years[1]:Years[2]
   N_dummyYears <- length(dummyYears)
   
   ## CREATE A DUMMY DATASET
   D_Id <- rep(dummyId, N_dummyYears)
   D_Year <- rep.int(dummyYears, rep(N_dummyId,N_dummyYears))
   dummy <- cbind.data.frame(Id = D_Id, Year = D_Year)
   
   ## CREATE THE COMBINED DETECTION ARRAY OF ALIVE & DUMMY DETECTIONS
   dummy <- rbind(dummy, data@data[ ,c("Id", "Year")])
   tab <- table(dummy$Id, dummy$Year)
   tab.ar <- matrix(0, dim(tab)[1], dim(tab)[2])
   tab.ar[ ,2:dim(tab)[2]] <- tab[ ,1:(dim(tab)[2]-1)] - 1
   tab.ar[tab.ar > 0] <- 1

   dimnames(tab.ar) <- dimnames(tab)

   return(tab.ar)
}

