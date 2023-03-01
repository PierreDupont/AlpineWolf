rm(list=ls())



library(Hmisc)
library(corrplot)
## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")

## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.5.2"
thisDir <- file.path(analysisDir, modelName)



flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


load(file.path(thisDir,"Habitat.RData"))
load(file.path(thisDir,"Detectors.RData"))

# names(habitat$grid)
det_cov <- detectors$grid[,c(3,5,7,16,20)]
det_cov <- st_drop_geometry(det_cov)
colnames(det_cov) <- c("Transects_L","Transects_Exp", "Snowfall", "Human_Pop", "East/West")


hab_cov <- habitat$grid[,c(4,5,6,12,17)]
hab_cov <- st_drop_geometry(hab_cov)
colnames(hab_cov) <- c("Forest%","Herbaceous%", "BareRocks%", "Human_Pop", "WHP")


hab_mat<-rcorr(as.matrix(hab_cov))
flattenCorrMatrix(hab_mat$r, hab_mat$P)

det_mat<-rcorr(as.matrix(det_cov))
flattenCorrMatrix(det_mat$r, det_mat$P)


corrplot(hab_mat$r, type="upper", order="hclust", 
         p.mat = hab_mat$P, sig.level = 0.01, insig = "blank", diag=FALSE)

corrplot(det_mat$r, type="upper", order="hclust", 
         p.mat = det_mat$P, sig.level = 0.01, insig = "blank", diag=FALSE)


