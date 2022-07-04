###############################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------ CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------ IMPORT REQUIRED LIBRARIES ------
library(rgdal)
library(raster)
library(coda)
library(nimble)
library(nimbleSCR)
library(stringr)
library(abind)
library(R.utils)
library(sf)
library(fasterize)
library(dplyr)
library(lubridate)
library(stars)
library(RANN)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(gridExtra)
library(MetBrewer)


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.SubSample.ResultsOutput"
thisDir <- file.path(analysisDir, modelName)
sex <- c("female","male")
status <- c("alpha","pup","other")
HabCov <- c("Bare Rocks","Herbaceous", "Forest","Pop","Wolf presence")
DetCov <- c("Transects_L", "Transect_Exp", "Snow", "East/West", "Log Pop")

myCols <- matrix(met.brewer(name = "Isfahan1", n = 8, type = "discrete")[c(2:4,7:5)],
                 nrow = 3, ncol = 2, byrow = F)
myCols2 <- matrix(met.brewer(name = "Pissaro", n = 7, type = "discrete")[c(1:2,4,5:6)],
                 nrow = 5, ncol = 1, byrow = F)


## -----------------------------------------------------------------------------
## ------  SUB-SAMPLING FIGURES -----
scenarios <- c(25,50,75,100)
pdf(file = file.path(thisDir, "figure_parameters_subsamples.pdf"),
     width = 10, height = 3.5)

##----- PLOTS N -----
par(mfrow = c(1,2))

## MEAN
ylim <- 5000
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Posterior estimate",
     ylab = "Population size (N)", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0, ylim, length.out = 5),
     labels = seq(0, ylim,length.out = 5))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.", scenarios[sc],"_mcmc.RData")))
  
  plot.violins2( dat.list = list(res$sims.list$N),
                 x = sc,
                 at = sc,
                 add = T,
                 col = "dodgerblue4",
                 violin.width = 0.05,
                 alpha = 0.9,
                 border.col = "dodgerblue4",
                 scale.width = F)
}#sc


  pdf(file = file.path(thisDir, "figure_parameters_subsamples.pdf"),
        width = 10, height = 3.5)
  
  
  par(mfrow=c(1,2))
  ## MEAN
  ylim <- 5000

  plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Mean estimate",
       ylab = "Population size (N)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 5),
       labels = seq(0,ylim,length.out = 5))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
    plot.violins2( dat.list = list(res$sims.list$N),
                   x = sc,
                   at = sc,
                   add = T,
                   col = "dodgerblue4",
                   violin.width = 0.05,
                   alpha = 0.9,
                   border.col = "dodgerblue4",
                   scale.width = F)
  }#sc
  
  plot(1,1, xlim = c(0,5), ylim = c(0,0.5),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = "CV(N)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,0.5,length.out = 6),
       labels = seq(0,0.5,length.out = 6))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
    points( y = res$sd$N/res$mean$N,
            x = sc,
            pch = 21, cex = 2,
            bg = "dodgerblue4",
            col = "dodgerblue4")
  }#sc


## CV
ylim <- 0.3
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Coefficient of variation",
     ylab = "CV(N)", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0,ylim,length.out = 5),
     labels = seq(0,ylim,length.out = 5))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
  points( y = res$sd$N/res$mean$N,
          x = sc,
          pch = 21, cex = 2,
          bg = "dodgerblue4",
          col = "dodgerblue4")
}#sc

  

##----- PLOTS sigma -----
par(mfrow=c(1,2))

## MEAN
ylim <- 20
valNoiseStatus <- c(-0.25,0,0.25)
valNoiseSex <- c(-0.05,0.05)
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Posterior estimate",
     ylab = "Sigma", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0,ylim,length.out = 5),
     labels = seq(0,ylim,length.out = 5))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.", scenarios[sc], "_mcmc.RData")))
  for(s in 1:length(sex)){
    for(ss in 1:length(status)){
      temp <- res$sims.list$sigma[ ,ss,s]
      
      plot.violins2( dat.list = list(temp),
                     x = sc + valNoiseStatus[ss] + valNoiseSex[s],
                     at = sc + valNoiseStatus[ss] + valNoiseSex[s],
                     add = T,
                     col = myCols[ss,s],
                     violin.width = 0.05,
                     alpha = 1,
                     border.col = myCols[ss,s],
                     scale.width = F)
    }#status
  }#sex
}#scenarios
legend(title = "Females",
       legend = c("Alpha","Pup","Other"),
       x = 0, y = 8, fill = myCols[1:3], cex = 0.6)
legend(title = "Males",
       legend = c("Alpha","Pup","Other"),
       x = 1, y = 8, fill = myCols[4:6], cex = 0.6)


## CV
ylim <- 1
valNoiseStatus <- c(-0.25,0,0.25)
valNoiseSex <- c(-0.05,0.05)
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Coefficient of variation",
     ylab = "CV(sigma)", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0,ylim,length.out = 6),
     labels = seq(0,ylim,length.out = 6))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
  for(s in 1:length(sex)){
    for(ss in 1:length(status)){
          points( y = res$sd$sigma[ss,s]/res$mean$sigma[ss,s],
                  x = sc + valNoiseStatus[ss] + valNoiseSex[s],
                  pch = 21, cex = 0.9,
                  bg  = c(myCols[ss,s]),
                  col = c(myCols[ss,s]))
    }#status
  }#sex
}#scenarios



  
##----- PLOTS theta -----
par(mfrow=c(1,2))

## MEAN
ylim <- 1
valNoiseStatus <- c(-0.25,0,0.25)
valNoiseSex <- c(-0.05,0.05)
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Posterior estimate",
     ylab = "theta", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0,ylim,length.out = 5),
     labels = seq(0,ylim,length.out = 5))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.", scenarios[sc], "_mcmc.RData")))
  for(s in 1:length(sex)){
    for(ss in 1:length(status)){
      temp <- res$sims.list$theta[ ,ss,s]
      
      plot.violins2( dat.list = list(temp),
                     x = sc + valNoiseStatus[ss] + valNoiseSex[s],
                     at = sc + valNoiseStatus[ss] + valNoiseSex[s],
                     add = T,
                     col = myCols[ss,s],
                     violin.width = 0.05,
                     alpha = 1,
                     border.col = myCols[ss,s],
                     scale.width = F)
    }#status
  }#sex
}#scenarios
legend(title = "Females",
       legend = c("Alpha","Pup","Other"),
       x = 0, y = 8, fill = myCols[1:3], cex = 0.6)
legend(title = "Males",
       legend = c("Alpha","Pup","Other"),
       x = 1, y = 8, fill = myCols[4:6], cex = 0.6)


## CV
ylim <- 1
valNoiseStatus <- c(-0.25,0,0.25)
valNoiseSex <- c(-0.05,0.05)
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Coefficient of variation",
     ylab = "CV(theta)", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0,ylim,length.out = 6),
     labels = seq(0,ylim,length.out = 6))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
  for(s in 1:length(sex)){
    for(ss in 1:length(status)){
      points( y = res$sd$theta[ss,s]/res$mean$theta[ss,s],
              x = sc + valNoiseStatus[ss] + valNoiseSex[s],
              pch = 21, cex = 0.9,
              bg  = c(myCols[ss,s]),
              col = c(myCols[ss,s]))
    }#status
  }#sex
}#scenarios



##----- PLOTS p0 -----
par(mfrow=c(1,2))

## MEAN
ylim <- 0.2
valNoiseStatus <- c(-0.25,0,0.25)
valNoiseSex <- c(-0.05,0.05)
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "p0 estimate",
     ylab = "p0", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0,ylim,length.out = 5),
     labels = seq(0,ylim,length.out = 5))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.", scenarios[sc], "_mcmc.RData")))
  for(s in 1:length(sex)){
    for(ss in 1:length(status)){
      temp <- res$sims.list$p0[ ,ss,s]
      
      plot.violins2( dat.list = list(temp),
                     x = sc + valNoiseStatus[ss] + valNoiseSex[s],
                     at = sc + valNoiseStatus[ss] + valNoiseSex[s],
                     add = T,
                     col = myCols[ss,s],
                     violin.width = 0.05,
                     alpha = 1,
                     border.col = myCols[ss,s],
                     scale.width = F)
    }#status
  }#sex
}#scenarios
legend(title = "Females",
       legend = c("Alpha","Pup","Other"),
       x = 0, y = 8, fill = myCols[1:3], cex = 0.6)
legend(title = "Males",
       legend = c("Alpha","Pup","Other"),
       x = 1, y = 8, fill = myCols[4:6], cex = 0.6)


## CV
ylim <- 2
valNoiseStatus <- c(-0.25,0,0.25)
valNoiseSex <- c(-0.05,0.05)
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Coefficient of variation",
     ylab = "CV(p0)", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0,ylim,length.out = 6),
     labels = seq(0,ylim,length.out = 6))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
  for(s in 1:length(sex)){
    for(ss in 1:length(status)){
      points( y = res$sd$p0[ss,s]/res$mean$p0[ss,s],
              x = sc + valNoiseStatus[ss] + valNoiseSex[s],
              pch = 21, cex = 0.9,
              bg  = c(myCols[ss,s]),
              col = c(myCols[ss,s]))
    }#status
  }#sex
}#scenarios





##----- PLOTS rho -----
par(mfrow = c(1,2))

## MEAN
ylim <- 1
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Posterior estimate",
     ylab = "% females (rho)", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0, ylim, length.out = 5),
     labels = seq(0, ylim,length.out = 5))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.", scenarios[sc],"_mcmc.RData")))
  
  plot.violins2( dat.list = list(res$sims.list$rho),
                 x = sc,
                 at = sc,
                 add = T,
                 col = "dodgerblue4",
                 violin.width = 0.05,
                 alpha = 0.9,
                 border.col = "dodgerblue4",
                 scale.width = F)
}#sc


## CV
ylim <- 0.3
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Coefficient of variation",
     ylab = "CV(rho)", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0,ylim,length.out = 7),
     labels = seq(0,ylim,length.out = 7))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
  points( y = res$sd$rho/res$mean$rho,
          x = sc,
          pch = 21, cex = 2,
          bg = "dodgerblue4",
          col = "dodgerblue4")
}#sc



##----- PLOTS psi -----
par(mfrow = c(1,2))

## MEAN
ylim <- 1
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Posterior estimate",
     ylab = "Inclusion probability (psi)", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0, ylim, length.out = 5),
     labels = seq(0, ylim,length.out = 5))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.", scenarios[sc],"_mcmc.RData")))
  
  plot.violins2( dat.list = list(res$sims.list$psi),
                 x = sc,
                 at = sc,
                 add = T,
                 col = "dodgerblue4",
                 violin.width = 0.05,
                 alpha = 0.9,
                 border.col = "dodgerblue4",
                 scale.width = F)
}#sc


## CV
ylim <- 0.3
plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
     type = "n", xaxt = "n", main = "Coefficient of variation",
     ylab = "CV(psi)", xlab = "", axes = F)
axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
axis(2, at = seq(0,ylim,length.out = 7),
     labels = seq(0,ylim,length.out = 7))
for(sc in 1:length(scenarios)){
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
  points( y = res$sd$psi/abs(res$mean$psi),
          x = sc,
          pch = 21, cex = 2,
          bg = "dodgerblue4",
          col = "dodgerblue4")
}#sc





## ---- PLOTS BetaDet -----
DetCov <- c("Transects_L", "Transect_Exp", "Snow", "East/West", "Log Pop")

for(b in 1:length(DetCov)){
  par(mfrow = c(1,2))
  
  ## MEAN
  ylim <- 2
  plot(1,1, xlim = c(0,5), ylim = c(-ylim,ylim),
       type = "n", xaxt = "n", main = "Posterior estimate",
       ylab = DetCov[b], xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(-ylim, ylim, length.out = 5),
       labels = seq(-ylim, ylim,length.out = 5))
  abline(h = 0, lty = 2)
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.", scenarios[sc],"_mcmc.RData")))
    
    plot.violins2( dat.list = list(res$sims.list$betaDet[,b]),
                   x = sc,
                   at = sc,
                   add = T,
                   col = "dodgerblue4",
                   violin.width = 0.05,
                   alpha = 0.9,
                   border.col = "dodgerblue4",
                   scale.width = F)
  }#sc
  
  ## CV
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[1],"_mcmc.RData")))
  ylim <- trunc(res$sd$betaDet[b]/abs(res$mean$betaDet[b]))+1
  plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = paste0("CV(",DetCov[b],")"), xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 7),
       labels = seq(0,ylim,length.out = 7))
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    points( y = res$sd$betaDet[b]/abs(res$mean$betaDet[b]),
            x = sc,
            pch = 21, cex = 2,
            bg = "dodgerblue4",
            col = "dodgerblue4")
  }#scenario
}#b



  
## ---- PLOTS BetaHab -----
HabCov <- c("Bare Rocks","Herbaceous", "Forest","Pop","Wolf presence")

for(b in 1:length(HabCov)){
  par(mfrow = c(1,2))
  
  ## MEAN
  ylim <- 3
  plot(1,1, xlim = c(0,5), ylim = c(-ylim,ylim),
       type = "n", xaxt = "n", main = "Posterior estimate",
       ylab = HabCov[b], xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(-ylim, ylim, length.out = 5),
       labels = seq(-ylim, ylim,length.out = 5))
  abline(h = 0, lty = 2)
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.", scenarios[sc],"_mcmc.RData")))
    
    plot.violins2( dat.list = list(res$sims.list$betaHab[,b]),
                   x = sc,
                   at = sc,
                   add = T,
                   col = "dodgerblue4",
                   violin.width = 0.05,
                   alpha = 0.9,
                   border.col = "dodgerblue4",
                   scale.width = F)
  }#sc
  
  ## CV
  load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[1],"_mcmc.RData")))
  ylim <- trunc(res$sd$betaHab[b]/abs(res$mean$betaHab[b]))+1
  plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = paste0("CV(",HabCov[b],")"), xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 7),
       labels = seq(0,ylim,length.out = 7))
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    points( y = res$sd$betaHab[b]/abs(res$mean$betaHab[b]),
            x = sc,
            pch = 21, cex = 2,
            bg = "dodgerblue4",
            col = "dodgerblue4")
  }#scenario
}#b


  
  plot(1,1, xlim = c(0,5), ylim = c(-3,ylim),
       type = "n", xaxt = "n", main = "Mean estimate",
       ylab = "BetaHab", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(-3,ylim,length.out = 5),
       labels = seq(-3,ylim,length.out = 5))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
    for(dc in 1:length(HabCov)){
      temp <- res$sims.list$betaHab[,dc]
    
    plot.violins2( dat.list = list(temp),
                   x = sc + valNoise[dc],
                   at = sc + valNoise[dc],
                   add = T,
                   col = myCols2[dc],
                   violin.width = 0.05,
                   alpha = 0.9,
                   border.col = myCols2[dc],
                   scale.width = F)
    }
  }#sc
  
  
  legend(title = "Habitat Covariates",
         legend = c("Bare Rocks","Herbaceous", "Forest","Human Pop","Wolf presence"),
         x=4, y=5, fill = myCols2[1:5], cex = 0.5)
  
#CV  
  
  plot(1,1, xlim = c(0,5), ylim = c(0,2.5),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = "CV(BetaHab)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,2.5,length.out = 6),
       labels = seq(0,2.5,length.out = 6))
  
  for(sc in 1:length(scenarios)) {
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
    for(dc1 in 1:length(DetCov)) {
      temp1 <- res$mean$betaHab[dc1]
      
      for(dc2 in 1:length(DetCov)) {
        temp2 <- res$sd$betaHab[dc2]
    points( y = temp2/temp1,
            x = sc,
            pch = 21, cex = 2,
            bg  = c(myCols2[dc1]),
            col = c(myCols2[dc1]))
  
      }
    }
  }#sc
  
  
   graphics.off()




## -----------------------------------------------------------------------------