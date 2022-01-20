#' @title Function to create input data for Efford's SECR from RovQuant input data. 
#' #'
#' @description
#' \code{SECRizeData} returns a list object with \code{} 
#' 
#' @param y Matrix of capture histories (rows: individuals; columns: detectors)
#' @param habitat.sp Spatial points dataframe with habitat grid points.
#' @param detector.sp Spatial points dataframe with detector locations.
#' @param plot.check Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated execution.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/SECRizeData.R
#' @keywords SECR
#'
#' @examples
#' # Create SECR input data from RovQuant data:
#' SECRizeData(...)
#' 
#' 
#'  



SECRizeData<-function(y,habitat.xy,detector.xy,detector.covs=NULL,detector.type="count",plot.check=TRUE){
## MELT DATA: CAPTURE HISTORY

y.melted <- melt(y, id=c())#--make sure y has proper dimnames


####  Keep ALL ROWS OF CAPTURES ONLY (detection > 0)

caps <- y.melted [which(y.melted $value>0),]


###### CREATE CAPTURE HISTORY DATA, MATCHING SECR FORMAT (CHECK THE EXACT COLUMNS, ORDER AND NAMES) #######

caps_df <- as.data.frame(caps)
colnames(caps_df) <- c("ID", "trapID", "freq")

## ADD TWO COLUMNS TO THE CAPTURE HISTORY FILE: SESSION AND OCCASION

Session <- 1984
Occasion <- 1

caps2 <- cbind(Session, caps_df, Occasion)



capfile <- caps2[,c("Session","ID","Occasion","trapID","freq")]
capfile<-capfile[rep(1:dim(capfile)[1],capfile$freq),]

capfile$freq <- NULL

####### CREATE trap-layout file

## CREATE A VECTOR OF TRAP ID AND ADD IT AS A COLUMN TO TRAPFILE


trapfile <- data.frame(trapID=1:dim(detector.xy)[1], detector.xy)


####### IMPORT TWO DATA FILES into SECR #######

temp <- merge(capfile,trapfile,all.x=TRUE)

table(temp$ID)

### CHANGE TO fmt 'XY' RATHER THAN 'trapID FORMAT, NOTE XY CASE SENSITIVITY

capt<-temp[,-1]
names(capt)<-c("Session" , "ID"   ,    "Occasion" ,"X"   ,     "Y" )



## ASSIGN OUR TRAP LAYOUT A NAME THAT SECR CAN USE AND PLOT IT

traps <- read.traps(data = trapfile, detector=detector.type)

##----- ADD COVARIATES TO TRAPS

if(!is.null(detector.covs))covariates(traps)<-data.frame(detector.covs)


## CHECK IF XY COORDINATE NUMBERS MATCH AND ARE OF THE EXACT SAME NUMBER OF DIGITS

#all(capt$X%in%traps$x)
#all(capt$Y%in%traps$y)

####### CREATE THE MODEL OBJECT: CAPTURE HISTORY AND PLOT IT  ########

capt <- make.capthist(capt, traps, fmt = "XY")#, noccasions = 1)#, covnames = NULL)


####### CREATE MASK OBJECT  ########

maskfile <- data.frame(maskID=1:dim(habitat.xy)[1], habitat.xy)

mask<-read.mask(data=maskfile)



#secr0 <- secr.fit(capt)
#secr0 <- secr.fit(capt,model=list(D~1,g0~1,sigma~1),mask=mask)
#secr0
##---IF PLOT

if(plot.check){
   plot(y~x,mask,pch=19,col="grey")
   points(traps$y~traps$x,col="navy",pch=19)
  }


##---- OUTPUT

out<-list(
   trap.secr=traps,
   mask.secr=mask,
   capt.secr=capt
   
)

return(out)

}



