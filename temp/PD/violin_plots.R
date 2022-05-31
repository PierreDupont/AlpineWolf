## ------  SUB-SAMPLING FIGURES -----
load(file.path(WD, "results.RData"))


scenarios <- c(0.25,0.5,0.75)

## PLOTS N
{ 
  graphics.off()
  pdf(file = file.path(WD, "figure_N.pdf"),
      width = 10, height = 3.5)

  ## MEAN
  ylim <- 3000
  
  plot(1,1, xlim = c(0,1.25), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Mean estimate",
       ylab = "Population size (N)", xlab = "", axes = F)
  axis(1, at = c(0.25,0.5,0.75,1), labels = c(0.25,0.5,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 5),
       labels = seq(0,ylim,length.out = 5))
  
  for(sc in scenarios){
    temp <- Results$N$mean[Results$N$scenario == sc, ]
    plot.violins2( dat.list = list(temp$mean),
                   x = sc,
                   at = sc,
                   add = T,
                   col = "dodgerblue4",
                   violin.width = 0.05,
                   alpha = 0.9,
                   border.col = "dodgerblue4",
                   scale.width = F)
  }#sc

  points( y = Results$N$mean[Results$N$scenario == 1.0, ],
          x = 1.0,
          pch = 21, cex = 2,
          bg = "dodgerblue4",
          col = "dodgerblue4")
  
  
  ## COEFFICIENT OF VARIATION
  ylim <- 1.0
  plot(1,1, xlim = c(0,1.25), ylim = c(0,ylim),
       type = "n", xaxt = "n", main = "Coefficient of variation",
       ylab = "CV(N)", xlab = "", axes = F)
  axis(1, at = c(0.25,0.5,0.75,1), labels = c(0.25,0.5,0.75,1))
  axis(2, at = seq(0,ylim,length.out = 5),
       labels = seq(0,ylim,length.out = 5))
  for(sc in scenarios){
    temp <- Results$N$CV[Results$N$scenario == sc, ]
    plot.violins2( dat.list = list(temp$mean),
                   x = sc,
                   at = sc,
                   add = T,
                   col = "dodgerblue4",
                   violin.width = 0.05,
                   alpha = 0.9,
                   border.col = "dodgerblue4",
                   scale.width = F)
  }#sc
  
  points( y = Results$N$CV[Results$N$scenario == 1.0, ],
          x = 1.0,
          pch = 21, cex = 2,
          bg = "dodgerblue4",
          col = "dodgerblue4")

  graphics.off()
}
