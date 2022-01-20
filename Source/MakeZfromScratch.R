MakeZfromScratch <- function(data.alive,
                             data.dead = NULL, 
                             samplingMonths = 1:12){
   
   ## MAKE SURE ONLY ONE DEAD RECOVERY PER INDIVIDUAL
   if(sum(duplicated(data.dead$Id))>0){stop("INDIVIDUALS CANNOT BE DEAD TWICE!!!!")}
   
   ## LIST ALL INDIVIDUALS
   data.alive$Id <- as.character(data.alive$Id)
   data.dead$Id <- as.character(data.dead$Id)
   dummyId <- unique(c(data.alive$Id, data.dead$Id))
   N_dummyId <- length(dummyId)
   
   ## LIST ALL YEARS
   Years <- range(unique(c(data.alive$Year, data.dead$Year)))
   dummyYears <- Years[1]:Years[2]
   N_dummyYears <- length(dummyYears)
   
   ## CREATE A DUMMY DATASET
   D_Id <- rep(dummyId, N_dummyYears*12)
   D_Year <- rep(rep.int(dummyYears, rep(N_dummyId,N_dummyYears)), 12)
   D_Months <- rep.int(1:12, rep(N_dummyId*N_dummyYears, 12))
   dummy <- cbind.data.frame(Id = D_Id, Year = D_Year, Month = D_Months)
   
   ## CREATE THE COMBINED DETECTION ARRAY OF ALIVE & DUMMY DETECTIONS
   dummyAlive <- rbind(dummy, data.alive@data[ ,c("Id","Year","Month")])
   tabAlive <- table(dummyAlive$Id, dummyAlive$Month, dummyAlive$Year)
   tabAlive.ar <- array(NA, c(dim(tabAlive)))
   tabAlive.ar[] <- tabAlive - 1
   tabAlive.ar[tabAlive.ar > 0] <- 1
   tabTotal <- tabAlive.ar
   
   if(!is.null(data.dead)){
      ## CREATE THE COMBINED DETECTION ARRAY OF DEAD & DUMMY RECOVERIES 
      dummyDead <- rbind(dummy, data.dead@data[ ,c("Id","Year","Month")])
      tabDead <- table(dummyDead$Id, dummyDead$Month, dummyDead$Year)
      tabDead.ar <- array(NA, c(dim(tabDead)))
      tabDead.ar[] <- tabDead - 1
      tabDead.ar[tabDead.ar > 0] <- 2
      tabTotal <- tabAlive.ar + tabDead.ar  
      }#if
  
   dimnames(tabTotal) <- dimnames(tabAlive)
   
   ## RE-ORDER THE MONTHS TO MATCH THE BIOLOGICAL YEAR CONSIDERED (= PUT SAMPLING MONTHS AT THE BEGINNING OF THE YEAR)
   allMonths <- 1:12
   outMonths <- allMonths[!(allMonths %in% samplingMonths)]
   tabTotal <- tabTotal[ ,c(samplingMonths,outMonths), ]
   
   ## RECONSTRUCT MONTHLY z BASED ON DETECTIONS & DEAD RECOVERIES
   z <- tabTotal
   for(i in 1:dim(z)[1]){
      # SUBSET DATA TO ID i
      temp <- z[i,,]
      # IF AT LEAST ONE DETECTION (SHOULD ALWAYS BE THE CASE)
      if(any(z[i,,] > 0)){
         # GET THE RANGE OF DETECTIONS (= WHEN ID i WAS ALIVE or RECOVERED)
         range.det <- range(which(z[i,,] > 0))
         # IF ID i IS RECOVERED DEAD 
         if(any(z[i,,] >= 2)){
               # GET MONTH OF RECOVERY
               death <- which(z[i,,] >= 2)
               # if(any(which(z[i,,] == 1) > death)){
               #    stop(paste("Individual", dimnames(z)[[1]][i], "was detected ALIVE AFTER being recovered DEAD!"))
               #    }
               # ID i WAS ALIVE BETWEEN THE FIRST DETECTION AND THE MONTH OF RECOVERY
               temp[range.det[1]:death] <- 1
               # IF RECOVERY HAPPENED BEFORE THE LAST CAPTURE OCCASION
               if(death < length(temp)){
                  # ID i WAS DEAD FROM THE CAPTURE OCCASION AFTER ITS DEAD RECOVERY ON
                  temp[(death+1):length(temp)] <- 2
               }
            # IF ID i WAS ONLY DETECTED ALIVE  
            }else{
               # ID i WAS ALIVE AT LEAST BETWEEN THE FIRST & LAST DETECTIONS
               temp[range.det[1]:range.det[2]] <- 1
            }
         }
      
      z[i,,] <- temp
   }
   z[z==0] <- NA
   dimnames(z) <- dimnames(tabTotal)  
   return(z)
}

