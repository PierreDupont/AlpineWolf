#' Density objects preparation
#'
#' R utility function to prepare input objects for the \code{getDensity} function (\code{\link{getDensity}}).
#'
#' The \code{getDensityInput} function is used in advance of density extraction with function \code{getDensity}.
#'
#' @param regions Raster denoting the region(s) for which density will be calculated, at the desired resolution.
#' @param habitat Raster of the habitat (as used in the SCR model). 
#' @param s Array of MCMC samples for individual activity center x- and y-coordinates 
#' (make sure these are in the same projection than the \code{habitat} and \code{regions} objects)
#' @param resolution resolution used to extract density, only used if \code{regions} is a Spatial Polygon.
#' @param plot.check A visualization option (if TRUE); 
#' displays the raster of the different regions that will be used to extract density,
#' as well as potential mismatches between the \code{habitat} and \code{regions} rasters.
#'
#' @return This function returns a list of objects:
#' \itemize{
#' \item regions.r : raster of the regions for which density will be calculated, with the desired resolution.
#' \item habitat.id: matrix of habitat grid cell IDs (NA for non-habitat; 0 for habitat not in the \code{regions} considered). 
#' \item regions.rgmx: matrix of dimensions number of regions * number of habitat cells containing
#'  binary indicators denoting to which region each habitat cell belongs.
#' \item sx: matrix of rescaled MCMC samples of x-coordinates of individual activity centers for use in the Rcpp \code{getDensity} function.  
#' \item sy: matrix of rescaled MCMC samples of y-coordinates of individual activity centers for use in the Rcpp \code{getDensity} function.
#' }
#'
#' @author Pierre Dupont
#'
#' @importFrom graphics plot points ## need to be checked once we agree on the function.
#'
#' @examples
#'
#' @export
getDensityInput <- function( regions, 
                             habitat,
                             s = NULL,
                             # resolution = NULL,
                             plot.check = TRUE
){
  ##-- 1. If a raster is provided as input, use its resolution.
  ##-- Else, turn spatial object 'regions' into a raster with the desired 
  ##-- output resolution.
  # if(class(regions) %in% c("RasterLayer", "RasterStack", "RasterBrick")){
  #   regions.r <- regions
  #   resolution <- res(regions.r)
  # } else {
  #   regions.r <- raster(extent(regions), resolution = resolution)
  #   regions.r <- rasterize(regions, regions.r)
  # }
  
  ##-- 1. Change the resolution of the habitat raster to the desired output resolution.
  if(sum(res(habitat) > res(regions))){
    habitat <- raster::disaggregate(habitat, fact = res(habitat)/res(regions))
  } else {
    habitat <- raster::aggregate(habitat, fact = res(habitat)/res(regions))
  }
  
  
  ##-- 2. Crop the habitat raster to the extent of the region of interest
  habitat <- crop(habitat, extent(regions))
  
  
  ##-- 3. Check if all the cells in the 'regions' raster are 
  ##-- considered as suitable habitat in the 'habitat' raster
  habitat[habitat[] == 0] <- NA
  # regions.r <- resample(regions.r, habitat, method = "ngb")
  if(any(!is.na(regions[]) & is.na(habitat[]))){
    warning("Some cells from the 'regions' raster are not considered as suitable habitat")
  }
  
  if(plot.check){
    probCells <- habitat
    probCells[] <- NA
    probCells[!is.na(regions[]) & is.na(habitat[])] <- 1
    plot(habitat, col = "gray60", legend = F) 
    plot(regions, add = T, legend = F)
    plot(probCells, add = T, col = adjustcolor("red", alpha.f = 0.7), legend = F)
  }
  
  
  ##-- 4. Combine 'habitat' and 'regions' rasters   
  ##-- in the current version, all regions cells not in the habitat are discarded.
  ##-- this means no density will be calculated for these cells later on.
  regions[is.na(regions[])] <- 0
  newRaster <- habitat + regions - 1
  
  
  ##-- 5. Create a matrix of regions binary indicators
  ## (rows == number of distinct regions ; columns == habitat raster cells)
  regionsNames <- sort(unique(na.omit(newRaster[])))
  regionsNames <- regionsNames[regionsNames > 0]
  regions.rgmx <- do.call(rbind, lapply(regionsNames, function(x)newRaster[] == x))
  regions.rgmx[is.na(regions.rgmx)] <- 0
  if(is.factor(regions)){##
    row.names(regions.rgmx) <- factorValues(regions, regionsNames)[,1]
    }else{
      row.names(regions.rgmx) <- regionsNames
    }
  
  ##-- Identify which cells are habitat
  isHabitat <- which(!is.na(newRaster[]))
  ##-- Keep only habitat cells  
  regions.rgmx <- regions.rgmx[ ,isHabitat]
  
  ##-- Force region.rgmx to be a matrix. 
  if(is.null(dim(regions.rgmx))){
    regions.rgmx <-  matrix(regions.rgmx, nrow = 1)
    rownames(regions.rgmx) <- regionsNames
  }
    
  ##-- 6. Create a matrix of habitat cell IDs
  cellIDs <- rep(NA, ncell(habitat))
  cellIDs[isHabitat] <- 1:length(isHabitat)
  habitat.id <- matrix( data = cellIDs,
                        nrow = dim(regions)[1],
                        ncol = dim(regions)[2],
                        byrow = TRUE)
  
  ##-- 7. Retrieve habitat x and y center
  
  
  ##-- 8. Rescale activity center x- and y- coordinates
  if(!is.null(s)){
    rescaled_s <- scaleCoordsToHabitatGrid( coordsData = s,
                                            coordsHabitatGridCenter = coordinates(regions),
                                            scaleToGrid = T)
    
    
    whereXY <- which(unlist(lapply(dimnames(rescaled_s$coordsDataScaled), function(x)any(x %in% c("x","y")))))
    split_s <- asplit(rescaled_s$coordsDataScaled, MARGIN = whereXY)
    
    output <- list( regions.r = newRaster,
                    habitat.id = habitat.id,
                    habitat.xy = rescaled_s$coordsHabitatGridCenterScaled[isHabitat,],
                    regions.rgmx = regions.rgmx,
                    sx = split_s$x,
                    sy = split_s$y)
    
  }else{
    rescaled_hab <- scaleCoordsToHabitatGrid( coordsData = coordinates(regions),
                                            coordsHabitatGridCenter = coordinates(regions),
                                            scaleToGrid = T)$coordsHabitatGridCenterScaled
    
    habitat.xy  <- rescaled_hab[isHabitat,]
    
    rescaled_s <- scaleCoordsToHabitatGrid( coordsData = s,
                                            coordsHabitatGridCenter = coordinates(regions),
                                            scaleToGrid = T)$coordsDataScaled
    
    output <- list( regions.r = newRaster,
                    habitat.id = habitat.id,
                    regions.rgmx = regions.rgmx,
                    habitat.xy = habitat.xy)
    
  }
  
  
  ##-- 9. OUTPUT LIST
  
  return(output)
}
