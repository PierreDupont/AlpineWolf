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
sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))


## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.SubSample.ResultsOutput"
thisDir <- file.path(analysisDir, modelName)
sex <- c("female","male")
status <- c("alpha","pup","other")
HabCov <- c("Bare Rocks","Herbaceous", "Forest","Pop","Wolf presence")
DetCov <- c("Transects_L", "Transect_Exp", "Snow", "East/West", "Log Pop")

myCols <- matrix(met.brewer(name="Isfahan1",n=8,type="discrete")[c(2:4,7:5)],
                 nrow = 3, ncol=2, byrow = F)
myCols2 <- matrix(met.brewer(name="Pissaro",n=7,type="discrete")[c(1:2,4,5:6)],
                 nrow = 5, ncol=1, byrow = F)

## ------  SUB-SAMPLING FIGURES -----


scenarios <- c(25,50,75,100)

##----- PLOTS N -----


  # pdf(file = file.path(thisDir, "figure_parameters_subsamples.pdf"),
  #     width = 10, height = 3.5)
  
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


  

  
  # ## COEFFICIENT OF VARIATION
  # ylim <- 1.0
  # plot(1,1, xlim = c(0,1.25), ylim = c(0,ylim),
  #      type = "n", xaxt = "n", main = "Coefficient of variation",
  #      ylab = "CV(N)", xlab = "", axes = F)
  # axis(1, at = c(25,50,75,100), labels = c(25,50,75,100))
  # axis(2, at = seq(0,ylim,length.out = 5),
  #      labels = seq(0,ylim,length.out = 5))
  # for(sc in scenarios){
  #   temp <- Results$N$CV[Results$N$scenario == sc, ]
  #   plot.violins2( dat.list = list(temp$mean),
  #                  x = sc,
  #                  at = sc,
  #                  add = T,
  #                  col = "dodgerblue4",
  #                  violin.width = 0.05,
  #                  alpha = 0.9,
  #                  border.col = "dodgerblue4",
  #                  scale.width = F)
  # }#sc
  # 
  # points( y = Results$N$CV[Results$N$scenario == 1.0, ],
  #         x = 1.0,
  #         pch = 21, cex = 2,
  #         bg = "dodgerblue4",
  #         col = "dodgerblue4")
  

##----- PLOTS sigma -----

  
  par(mfrow=c(1,2))
  ## MEAN
  
 ylim <- 8
 valNoise <- c(-0.2,0,0.2)

 
 plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Sigma estimate",
       ylab = "Sigma", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 5),
       labels = seq(0,ylim,length.out = 5))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
 
    for(s in 1:length(sex)){
      for(ss in 1:length(status)){
        temp <- res$sims.list$sigma[ ,ss,s]
       
    plot.violins2( dat.list = list(temp),
                   x = sc + valNoise[ss],
                   at = sc + valNoise[ss],
                   add = T,
                   col = myCols[ss,s],
                   violin.width = 0.05,
                   alpha = 1,
                   border.col = myCols[ss,s],
                   scale.width = F)
    }#status
  }#sex
}
  legend(title = "Females",
          legend = c("Alpha","Pup","Other"),
         x = 0, y = 8, fill = myCols[1:3], cex = 0.6)
   legend(title = "Males",
          legend = c("Alpha","Pup","Other"),
          x = 1, y = 8, fill = myCols[4:6], cex = 0.6)
  
  
  
  plot(1,1, xlim = c(0,5), ylim = c(0,1),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = "CV(sigma)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,1,length.out = 6),
       labels = seq(0,1,length.out = 6))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
 
    for(s in 1:length(sex)){
      for(ss in 1:length(status)){
        temp1 <- res$mean$sigma[ss,s]
        
    for(q in 1:length(sex)){
          for(qq in 1:length(status)){
            temp2 <- res$sd$sigma[qq,q]

    points( y = temp2/temp1,
            x = sc + valNoise[ss],
            pch = 21, cex = 2,
            bg  = c(myCols[ss,s]),
            col = c(myCols[ss,s]))
        }
       }
      }
     }
    }
  
  
  
  # ## COEFFICIENT OF VARIATION
  # ylim <- 1.0
  # plot(1,1, xlim = c(0,1.25), ylim = c(0,ylim),
  #      type = "n", xaxt = "n", main = "Coefficient of variation",
  #      ylab = "CV(N)", xlab = "", axes = F)
  # axis(1, at = c(25,50,75,100), labels = c(25,50,75,100))
  # axis(2, at = seq(0,ylim,length.out = 5),
  #      labels = seq(0,ylim,length.out = 5))
  # for(sc in scenarios){
  #   temp <- Results$N$CV[Results$N$scenario == sc, ]
  #   plot.violins2( dat.list = list(temp$mean),
  #                  x = sc,
  #                  at = sc,
  #                  add = T,
  #                  col = "dodgerblue4",
  #                  violin.width = 0.05,
  #                  alpha = 0.9,
  #                  border.col = "dodgerblue4",
  #                  scale.width = F)
  # }#sc
  # 
  # points( y = Results$N$CV[Results$N$scenario == 1.0, ],
  #         x = 1.0,
  #         pch = 21, cex = 2,
  #         bg = "dodgerblue4",
  #         col = "dodgerblue4")
  
  
  
##----- PLOTS theta -----
  
  
  ylim <- 1
  par(mfrow=c(1,2))
  ## MEAN

  
  plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Mean estimate",
       ylab = "Theta", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 5),
       labels = seq(0,ylim,length.out = 5))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
    for(s in 1:length(sex)){
      for(ss in 1:length(status)){
        temp <- res$sims.list$theta[ ,ss,s]
        
        plot.violins2( dat.list = list(temp),
                       x = sc + valNoise[ss],
                       at = sc + valNoise[ss],
                       add = T,
                       col = myCols[ss,s],
                       violin.width = 0.05,
                       alpha = 1,
                       border.col = myCols[ss,s],
                       scale.width = F)
      }#status
    }#sex 
  } #sc
  legend(title = "Females",
         legend = c("Alpha","Pup","Other"),
         x = 0.5, y = 1, fill = myCols[1:3], cex = 0.5)
  legend(title = "Males",
         legend = c("Alpha","Pup","Other"),
         x = 0.8, y = 1, fill = myCols[4:6], cex = 0.5)
  
  
  
  plot(1,1, xlim = c(0,5), ylim = c(0,1),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = "CV(theta)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,1,length.out = 6),
       labels = seq(0,1,length.out = 6))
  
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
  for(s in 1:length(sex)){
    for(ss in 1:length(status)){
      temp1 <- res$mean$theta[ ss,s]
      
      for(q in 1:length(sex)){
        for(qq in 1:length(status)){
          temp2 <- res$sd$theta[ qq,q]
          
        
    points( y = temp2/temp1,
            x = sc + valNoise[ss,s],
            pch = 21, cex = 2,
            bg = c(myCols[ss,s], myCols[qq,q]),
            col = c(myCols[ss,s], myCols[qq,q]))
 
        }
      }
    }
  }
}#sc
  
  
## ---- PLOTS p0----
  
  par(mfrow=c(1,2))
  
  ylim <- 0.15
  ## MEAN
  
  plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Mean estimate",
       ylab = "p0", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 5),
       labels = seq(0,ylim,length.out = 5))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
    for(s in 1:length(sex)){
      for(ss in 1:length(status)){
        temp <- res$sims.list$p0[ ,ss,s]
        
        plot.violins2( dat.list = list(temp),
                       x = sc + valNoise[ss],
                       at = sc + valNoise[ss],
                       add = T,
                       col = myCols[ss,s],
                       violin.width = 0.05,
                       alpha = 1,
                       border.col = myCols[ss,s],
                       scale.width = F)
      }#status
    }#sex 
  } #sc
  
  legend(title = "Females",
         legend = c("Alpha","Pup","Other"),
         x=3.5, y=0.15, fill = myCols[1:3], cex = 0.5)
  legend(title = "Males",
         legend = c("Alpha","Pup","Other"),
         x=4.3, y=0.15, fill = myCols[4:6], cex = 0.5)
  
  #CV
  
  plot(1,1, xlim = c(0,5), ylim = c(0,1.5),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = "CV(p0)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,1.5,length.out = 6),
       labels = seq(0,1.5,length.out = 6))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
    for(s in 1:length(sex)){
      for(ss in 1:length(status)){
        temp1 <- res$mean$p0[ ss,s]
        
        for(q in 1:length(sex)){
          for(qq in 1:length(status)){
            temp2 <- res$sd$p0[ qq,q]
            
            
            points( y = temp2/temp1,
                    x = sc,
                    pch = 21, cex = 2,
                    bg = c(myCols[ss,s], myCols[qq,q]),
                    col = c(myCols[ss,s], myCols[qq,q]))
            
          }
        }
      }
    }
  }#sc
  
## ---- PLOTS rho----
  
  
  par(mfrow=c(1,2))
  
  ylim <- 1
  ## MEAN

  plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Mean estimate",
       ylab = "rho", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 5),
       labels = seq(0,ylim,length.out = 5))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
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
  
  plot(1,1, xlim = c(0,5), ylim = c(0,0.5),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = "CV(rho)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,0.5,length.out = 6),
       labels = seq(0,0.5,length.out = 6))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))

        points( y = res$sd$rho/res$mean$rho,
            x = sc,
            pch = 21, cex = 2,
            bg = "dodgerblue4",
            col = "dodgerblue4")
  }#sc
  
## ---- PLOTS psi-----
  
  ylim <- 1
  par(mfrow=c(1,2))
  ## MEAN

  
  plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Mean estimate",
       ylab = "psi", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 5),
       labels = seq(0,ylim,length.out = 5))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
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
  
  plot(1,1, xlim = c(0,5), ylim = c(0,0.5),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = "CV(psi)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,0.5,length.out = 6),
       labels = seq(0,0.5,length.out = 6))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))

        points( y = res$sd$psi/res$mean$psi,
            x = sc,
            pch = 21, cex = 2,
            bg = "dodgerblue4",
            col = "dodgerblue4")
  }#sc
  
  
## ---- PLOTS BetaDet -----
  par(mfrow=c(1,2))
  
  ylim <- 5
  ## MEAN
  valNoise <- c(-0.05,0,0.05)
  
  plot(1,1, xlim = c(0,4), ylim = c(-4,ylim),
       type = "n", xaxt = "n", main = "Mean estimate",
       ylab = "BetaDet", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(-4,ylim,length.out = 5),
       labels = seq(-4,ylim,length.out = 5))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
   
     for(dc in 1:length(DetCov)){
      temp <- res$sims.list$betaDet[,dc]
    
    plot.violins2( dat.list = list(temp),
                   x = sc + valNoise[,dc],
                   at = sc + valNoise[,dc],
                   add = T,
                   col = "dodgerblue4",
                   violin.width = 0.02,
                   alpha = 0.9,
                   border.col = "dodgerblue4",
                   scale.width = F)
    }#habcov
  }#sc
  
  plot(1,1, xlim = c(0,5), ylim = c(0,0.5),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = "CV(BetaDet)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,0.5,length.out = 6),
       labels = seq(0,0.5,length.out = 6))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    for(dc1 in 1:length(DetCov)){
        temp1 <- res$mean$betaDet[dc1]
        
        for(dc2 in 1:length(DetCov)){
            temp2 <- res$sd$betaDet[dc2]

        
    points( y = temp2/temp1,
            x = sc + valNoise[,dc],
            pch = 21, cex = 2,
            bg  = c(myCols2[,dc]),
            col = c(myCols2[,dc]))
        }
    }
}#sc
  
  
  
## ---- PLOTS BetaHab -----
  
  
  par(mfrow=c(1,2))
  ## MEAN

  
  plot(1,1, xlim = c(0,5), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Mean estimate",
       ylab = "BetaHab", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 5),
       labels = seq(0,ylim,length.out = 5))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
    plot.violins2( dat.list = list(res$sims.list$betaHab),
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
       ylab = "CV(BetaHab)", xlab = "", axes = F)
  axis(1, at = c(1,2,3,4), labels = c(0.25,0.50,0.75,1))
  axis(2, at = seq(0,0.5,length.out = 6),
       labels = seq(0,0.5,length.out = 6))
  
  for(sc in 1:length(scenarios)){
    load(file.path(thisDir, paste0("AlpineWolf.SubSample.",scenarios[sc],"_mcmc.RData")))
    
    for(i in seq_along(res$mean$betaHab))
      for(j in seq_along(res$sd$theta))
    points( y = res$sd$betaHab[j]/res$mean$betaHab[i],
            x = sc,
            pch = 21, cex = 2,
            bg  = c(myCols2[ss,s]),
            col = c(myCols2[ss,s])
  }#sc
  
  
  # graphics.off()


