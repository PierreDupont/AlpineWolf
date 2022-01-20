#' @title MakeAugmentation function
#'
#' @description \code{MakeAugmentation} increases the dimensions of an object along
#'  the individual and/or the year dimension. It returns a \code{Vector}, \code{Matrix} 
#'  or \code{Array} object with the expanded object.
#'
#' @param y A \code{Vector}, \code{Matrix} or \code{Array} object containing the individual detection histories.
#' @param aug.factor A \code{numeric} object defining the augmentation factor to be used.
#' @param replace.value A \code{numeric} object defining the value to be repeated for augmented individuals.
#' 
#' @return A \code{Vector}, \code{Matrix} or \code{Array} object containing the augmented y.

MakeAugmentation <- function( y,
                              aug.factor= NULL,
                              aug.years = NULL,
                              replace.value = NA){
  ## Vector Data augmentation
  if(is.vector(y)){
    if(is.null(names(y))){ 
      names(y) <- 1:length(y) 
    }
    
    if(!is.null(aug.factor)){
      y.aug <- c(y, rep(replace.value, round(length(y) * aug.factor)))
      names(y.aug) <- c(names(y), rep("Augmented", round(length(y) * aug.factor)))
      y <- y.aug
    }
    
    if(!is.null(aug.years)){
      y.aug <- c(y, rep(replace.value, aug.years))
      maxYear <- as.numeric(names(y)[length(names(y))])
      augmentedYears <- (maxYear + 1) : (maxYear + aug.years)
      names(y.aug) <- c(names(y), augmentedYears)
      y <- y.aug
    }
  }
  
  ## Matrix augmentation
  if(is.matrix(y)){
    if(is.null(dimnames(y))){
      dimnames(y) <- list(1:dim(y)[1], 1:dim(y)[2])
    }
    
    if(!is.null(aug.factor)){
      n.tot <- round(dim(y)[1]*(1 + aug.factor))
      y.aug <- matrix(replace.value, n.tot, dim(y)[2])
      y.aug[1:dim(y)[1], ] <- y
      dimnames(y.aug) <- list(c( dimnames(y)[[1]], rep("Augmented", n.tot - dim(y)[1])),
                              dimnames(y)[[2]])
      y <- y.aug
    }
    
    if(!is.null(aug.years)){
      y.aug <- matrix(replace.value, dim(y)[1], dim(y)[2] + aug.years)
      y.aug[ ,1:dim(y)[2]] <- y
      maxYear <- as.numeric(dimnames(y)[[2]][length(dimnames(y)[[2]])])
      augmentedYears <- (maxYear + 1) : (maxYear + aug.years)
      dimnames(y.aug) <- list( dimnames(y)[[1]],
                               c(dimnames(y)[[2]], augmentedYears))
      y <- y.aug
    }
  }  
  
  ## 3D Array augmentation
  if(length(dim(y)) == 3){
    if(is.null(dimnames(y))){
      dimnames(y) <- list(1:dim(y)[1], 1:dim(y)[2], 1:dim(y)[3])
    }
    
    if(!is.null(aug.factor)){
      n.tot <- round(dim(y)[1]*(1 + aug.factor))
      y.aug <- array(replace.value, c(n.tot,dim(y)[2],dim(y)[3]))
      y.aug[1:dim(y)[1], , ] <- y
      dimnames(y.aug) <- list(c( dimnames(y)[[1]],rep("Augmented", n.tot - dim(y)[1])),
                              dimnames(y)[[2]],
                              dimnames(y)[[3]])
      y <- y.aug
    }
    
    if(!is.null(aug.years)){
      y.aug <- array(replace.value, c(dim(y)[1],dim(y)[2],dim(y)[3] + aug.years))
      y.aug[ , ,1:dim(y)[3]] <- y
      maxYear <- as.numeric(dimnames(y)[[3]][length(dimnames(y)[[3]])])
      augmentedYears <- (maxYear + 1) : (maxYear + aug.years)
      dimnames(y.aug) <- list( dimnames(y)[[1]],
                               dimnames(y)[[2]],
                               c(dimnames(y)[[3]], augmentedYears))
      y <- y.aug
    }
  }  
  
  ## 4D Array augmentation
  if(length(dim(y)) == 4){
    if(is.null(dimnames(y))){
      dimnames(y) <- list(1:dim(y)[1], 1:dim(y)[2], 1:dim(y)[3], 1:dim(y)[4])
    }
    
    if(!is.null(aug.factor)){
      n.tot <- round(dim(y)[3]*(1 + aug.factor))
      y.aug <- array(replace.value, c(dim(y)[1], dim(y)[2], n.tot, dim(y)[4]))
      y.aug[ , ,1:dim(y)[3], ] <- y
      dimnames(y.aug) <- list( dimnames(y)[[1]],
                               dimnames(y)[[2]],
                               c(dimnames(y)[[3]], rep("Augmented",  n.tot - dim(y)[3])),
                               dimnames(y)[[4]])
      y <- y.aug
    }
    
    if(!is.null(aug.years)){
      y.aug <- array(replace.value, c(dim(y)[1], dim(y)[2], dim(y)[3], dim(y)[4] + aug.years))
      y.aug[ , , ,1:dim(y)[4]] <- y
      maxYear <- as.numeric(dimnames(y)[[4]][length(dimnames(y)[[4]])])
      augmentedYears <- (maxYear + 1) : (maxYear + aug.years)
      dimnames(y.aug) <- list( dimnames(y)[[1]],
                               dimnames(y)[[2]],
                               dimnames(y)[[3]],
                               c(dimnames(y)[[4]], augmentedYears))
      y <- y.aug
    }
  } 
  
  return (y)
}