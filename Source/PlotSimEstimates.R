PlotSimEstimates <- function( y1
                            , y2 = NULL   
                            , x 
                            , Levels = NULL
                            , simul.value = 0
                            , bias = FALSE
                            , relative = FALSE
                            , x.text = ""
                            , y.text = ""
                            , y.lab = NULL
                            , at.y = NULL
                            , x.lab = NULL
                            , at.x = NULL
                            , myColors = NULL
                            , min.TradeOff = FALSE)
   {
   y <- y1
   if(relative & !bias){y <- y1/simul.value}
   if(bias & !relative){y <- y1-simul.value}
   if(bias & relative){y <- (y1-simul.value)/simul.value}
   if(!is.null(y2)){y <- (y1-simul.value)/simul.value + y2/simul.value}
   
   ## Generate the different levels combinations
   if(is.null(Levels)){Levels <- list(rep(1,length(y)))}
   L <- lapply(Levels, unique)
   parm <- expand.grid(L)  
   parm <- as.matrix(parm)
   colnames(parm) <- names(Levels)
   
   ## Convert list of variables to a dataframe
   df <- do.call(cbind.data.frame, Levels)
   colnames(df) <- names(Levels)
   
   ## list of changing plot breakpoints
   change.plot <- seq(1, dim(parm)[1], length(L[[1]]))
   
   ## Set colors
   if(is.null(myColors)){myColors <- rep("deepskyblue4",length(L[[1]]))}
   myCols1 <- rep(myColors, length(change.plot))
   myCols <- as.data.frame(t(col2rgb(myCols1))/255)
   myCols$Color <- myCols1
   this.col <- NULL
   
   ## plot the panels
   for(l in 1:dim(parm)[1]){
      index <- which(unlist(apply(df,1, function(x)all(x == parm[l,]))))
   
      temp.x <- x[index]
      temp.y <- y[index]
      mean.y <- tapply(temp.y, temp.x, function(x){mean(x, na.rm = TRUE)})
      sd.y <- tapply(temp.y, temp.x,function(x){sd(x, na.rm = TRUE)})
      sd.y[is.na(sd.y)] <- 0
      lci <- mean.y - sd.y
      uci <- mean.y + sd.y
      
      if(l %in% change.plot){
         ## Set x & y labels
         if(is.null(x.lab)){x.lab <- as.numeric(as.character(names(mean.y)))}
         if(is.null(at.x)){at.x <- as.numeric(as.character(names(mean.y)))}
         if(is.null(y.lab)){
            ymin <- min(y)
            ymax <- max(y)
            y.lab <- round(seq(ymin, ymax,length.out = 10))
            }

         plot(as.numeric(as.character(names(mean.y))), mean.y, type = "n", axes = FALSE, xlab = "", ylab = "", ylim = c(min(at.y), max(at.y)), xlim = range(at.x))
         axis(1, at = at.x, labels = x.lab, tck = 0.01, cex.axis = 1.2, las = 1, hadj = 0.5)
         axis(2, at = y.lab, labels = y.lab, cex.axis = 1.2, las = 1, tck = 0.01)
         mtext(x.text, 1, 2.2, cex = 1.2, font = 2)
         mtext(y.text, 2, 2.6, cex = 1.2, font = 2)
         if(bias){abline(h = 0,lty = 2)}
         }
      
      polygon(c(as.numeric(as.character(names(mean.y[!is.na(mean.y)]))), rev(as.numeric(as.character(names(mean.y[!is.na(mean.y)]))))), 
              c(lci[!is.na(mean.y)], rev(uci[!is.na(mean.y)])), border = F, col = rgb(red = myCols[l,1], green = myCols[l,2], blue = myCols[l,3], alpha = 0.3))
      
      points(mean.y ~ as.numeric(as.character(names(mean.y))), type="l", lwd=2, col=rgb(red = myCols[l,1], green = myCols[l,2], blue = myCols[l,3], alpha = 1))
      
      this.col <- c(this.col, myCols[l,4])
      
      if(min.TradeOff){
         myMin <- which(mean.y == min(mean.y))
         abline(v = as.character(names(mean.y))[myMin],lty = 2, lwd=2, col = rgb(red = myCols[l,1], green = myCols[l,2], blue = myCols[l,3], alpha = 1))}
         }
   
   L.num <- do.call(paste,lapply(Levels, function(x) x))
   L.num2 <- unique(L.num)
   parm <- as.data.frame(parm)
   parm$n.rep <- unlist(lapply(L.num2, function(x){cbind(length(which(L.num==x)))}))
   parm$this.col <- this.col
   print(parm)
   }
