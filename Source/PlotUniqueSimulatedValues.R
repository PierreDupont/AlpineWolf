
#' @title Plot estimated values and 95%CI (x-axis) for each unique simulation repetitions (y-axis)
#'
#' @description
#' \code{PlotUniqueSimualedValues} GPlot estimated values and 95%CI (x-axis) for each unique simulation set
#' 
#' @param point.estimates A \code{vector} with a parameter point estimates (mean, median) obtained for each simulations set
#' @param CI.low A \code{vector} with a parameter lower Confidence intervals for each simulation set
#' @param CI.high A \code{vector} with a parameter higher Confidence intervals for each simulation set
#' @param simulated.values A \code{vector} with the value used to simulate the data set for each simulation set
#' @param ylim x limits
#' @param xlim y limits
#' @param main title of the graphic
#' @return A \code{list} A plot
#' @example 
#' 
#'                PlotUniqueSimualedValues( point.estimates = point.estimates
#'                                        , CI.low = CI.low
#'                                        , CI.high = CI.high
#'                                        , simulated.values = simulated.values
#'                                        , ylim = NULL
#'                                        , xlim = NULL
#'                                        , main = NULL
#'                                        )
#' 

PlotUniqueSimulatedValues <- function(  point.estimates = point.estimates
                                     , CI.low = CI.low
                                     , CI.high = CI.high
                                     , simulated.values = simulated.values
                                     , ylim = NULL
                                     , xlim = NULL
                                     , main = NULL 
                                     , col = NULL
                                     , line.col = NULL){

### assess the x and ylim 
if(is.null(xlim)){   
xlim=c(min(c(CI.low,CI.high)), max(c(CI.low,CI.high)))
}

nb.simul <-length(point.estimates)
if(is.null(ylim)){   
   ylim=c(0, nb.simul)
}

# plot empty graph
plot(1, type="n", xlab="", ylab="Simulations",
     xlim=xlim,
     ylim=ylim,
     main = main)


# get unique colors for different simulated values 
u.values <- unique(simulated.values)
nb.values <- length(u.values)
col.values <- as.numeric(as.factor(simulated.values))

if(is.null(col)){
   col <- c("red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet","red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet","red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet","red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet","red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet" )  
   col <- col[col.values]
}

if(is.null(line.col)){
   line.col <- c("red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet","red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet","red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet","red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet","red", "black", "blue", "green", "purple", "pink","yellow", "orange",
            "brown","gray29","red1" ,"azure3","springgreen1","darkviolet" )  
   line.col <- line.col[col.values]
}

## plot point estimates 
points((1:nb.simul) ~ point.estimates,pch=16, col=col)
## plot Confidence intervals 
segments(x0= CI.low 
         ,x1= CI.high
         ,y0= c(1:nb.simul)
         ,y1= c(1:nb.simul)
         ,col=col)

### add a vertical line for the true value
segments(x0= simulated.values 
         ,x1= simulated.values
         ,y0= c(1:nb.simul)-0.5
         ,y1= c(1:nb.simul)+0.5
         ,col=line.col)
}


