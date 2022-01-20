#' @title Function to convert a raster covariate to a vector that follows the dinhomPP requirements
#'
#' @description
#' \code{RasterToDinhomPP} returns a matrix with the columns corresponding to the different rasters. 
#' A column (a vector) can be used as an input for the covariates on density with the dinhomPP in Nimble. 
#' @param raster A \code{RasterLayer} object or rasterstack object when several covariates are to be used.
#'
#' @examples
#' # Generate a spatial point data frame with simulated activity center locations
#' dinhomPPcovs <- RasterToDinhomPP(raster)

RasterToDinhomPP <- function(raster){
   
   habQualityDims <- dim(raster)[2:1]
   HabQualityCovariate <- matrix(NA,prod(habQualityDims),dim(raster)[3])
   habQualityCovariate.arr <- array(NA, dim = dim(raster) )
   for( i in 1:dim(raster)[3]){
      # Convert raster to a matrix 
      habQualityCovariate.arr[,,i] <- as.matrix(raster[[i]])
      #habQualityCovariate.arr[,,i] <- habQualityCovariate.arr[dim(habQualityCovariate.arr)[1]:1,,i]
      HabQualityCovariate[,i] <- intensityRearrange(habQualityCovariate.arr[,,i], habQualityDims, c(2, 1), c(0, 0))
   } 
   
   toggle.habitat <- as.numeric(!is.na(HabQualityCovariate[,1]))
   
   HabQualityCovariate[is.na(HabQualityCovariate)] <- 0
   colnames(HabQualityCovariate) <- names(raster)
   return(list(HabQualityCovariate = HabQualityCovariate,
               toggle.habitat = toggle.habitat
   ))
   
}


