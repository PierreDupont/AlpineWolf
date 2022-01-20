#' @title Function to create a least-cost path that joins two points across a topography.
#'
#' @description
#' \code{SimulateDetection} returns a list object with \code{} 
#' 
#' @param p0 Numeric variable denoting the intercept value of the negative binomial detection function describing the decay of detection with increasing distance from the AC.
#' @param sigma Numeric variable denoting the scale parameter of the detection function.
#' @param sd.eps.p Numeric variable denoting the hyperparameter associated with the spatial association of detection locations for members of the same herd.
#' @param AC.sp Spatial points dataframe with herd activity center locations.
#' @param detector.sp Spatial points dataframe with detector locations.
#' @param covs A data frame with the covariate values associated with each habitat grid cell (order identical to "covs").
#' @param betas A numeric vector denoting the coefficient values associated with the covariates to be used during simulations (effect on detection probability).
#' @param plot Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated during simulations.
#' areas on.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/SimulateHerdDetection.R
#' @keywords simul
#'
#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' det.sim<-SimulateHerdDetection(p0=0.6,sigma=1200,AC.sp=AC.sp,detector.sp=detector.sp,covs=covs, betas=betas,plot=TRUE)
#' 
#'  



MakeLCPath<-function(x,A,B,plot=TRUE){
   
   AtoB.sl<-NULL
   AtoB.cost<-0
   A<-unlist(A)
   B<-unlist(B)
   
   if(plot){
      points(A[1],A[2],col="red",pch=19,cex=1.3)
      points(B[1],B[2],col="white",pch=19,cex=1.3)
      
   }
   
   try(
      {
         AtoB.sl <- shortestPath(x, A , B, output="SpatialLines")
         if(plot)lines(AtoB.sl, col="red", lwd=2,lend=1)
         AtoB.cost<-costDistance(x, A, B)
      },
      silent=TRUE
      )
   out<-list(AtoB.cost=AtoB.cost,
             AtoB.sl=AtoB.sl)

   }