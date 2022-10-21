getDensity <- nimbleFunction(run = function(
    s = double(2),
    habitatGrid = double(2),
    indicator = double(1),
    numWindows = double(0),
    n.individuals = double(0)){
  returnType(double(1))
  dens <- numeric(length = numWindows, value = 0)
  for(i in 1:n.individuals){
    if(indicator[i]==1){
      windowInd <- habitatGrid[trunc(s[i,2]) + 1, trunc(s[i,1]) + 1]
      dens[windowInd] <- 1 + dens[windowInd]
    }
  }
  return(dens)
})
