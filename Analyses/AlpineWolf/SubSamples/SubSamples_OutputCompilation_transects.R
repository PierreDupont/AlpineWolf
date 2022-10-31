###############################################################################
##### ----------------------- ALPINE WOLF SCR ---------------------------- #####
##### ----------------- s[CLC + HumPop + Zone] --------------------------- #####
##### ---------------- z[psi,rho,theta] ---------------------------------- #####
##### ---------------- y[p0(transects + zone + ),sigma] ------------------ #####
################################################################################
## ------1. CLEAN THE WORK ENVIRONMENT ------
rm(list=ls())


## ------2. IMPORT REQUIRED LIBRARIES ------
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
library(fs)
library(purrr)
library(ggplot2)
library(reshape2)
library(ghibli)
library(wesanderson)
library(hrbrthemes)


## ------ 3. SET REQUIRED WORKING DIRECTORIES ------
source("workingDirectories.R")


## ------ 4.  SOURCE CUSTOM FUNCTIONS ------
sourceDirectory(path = file.path(getwd(),"Source"), modifiedOnly = F)
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetDensity.cpp"))
#sourceCpp(file = file.path(getwd(),"Source/cpp/GetSpaceUse.cpp"))
##function to extract means ecc..
col_mean <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- mean(x[[i]])
  }
  output
}
col_sd <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- sd(x[[i]])
  }
  output
}
col_cv <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- sd(x[[i]]/mean(x[[i]]))
  }
  output
}
col_lci <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- qs(x[[i]],0.025)
  }
  output
}
col_uci <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- qs(x[[i]],0.975)
  }
  output
}
qs <- function(x,y){as.numeric(quantile(x,y))}


## -----------------------------------------------------------------------------
## ------ 5. SET ANALYSIS CHARACTERISTICS -----
## MODEL NAME 
modelName = "AlpineWolf.SubSample.ResultsOutput"
thisDir <- file.path(analysisDir, modelName)



## -----------------------------------------------------------------------------
## ------   6. PROCESS OUTPUTS -----

OutDir <- file.path(thisDir, "output")




resLista <-list()  #create a list which will contain all your results


rpp <- 100
sim_names <- c("25","50","75")  #Here the names of your simulation  (i.e. 25,50,75,100)

for(sc in 1:length(sim_names)){
  
  tempLista <-list() #temporary list for each scenario
  
  for(rp in 1:rpp){  #This will be from 1:50
    
    # for(rp in 6:rpp){  #This will be from 1:50
    # tryCatch({ 
    ## List output files for scenario "sc" and repetition "rp" only
    outputs <- list.files(OutDir, paste0("AlpineWolf.SubSample.Transects_",sim_names[sc],"_",rp,"_"))
    print(outputs)
    # if (file.size(outputs) < 0) stop("file is 0 KB!")
    # }, error= function(e){cat("ERROR:", conditionMessage(e), "\n")})
    ## Maybe add a check to make sure all files listed in "outputs" are bigger than 0 KB (but I don't know how to do this?)
    
    ## Collect bites from the different chains for this percentage and this repetition only
    nimOutput <- collectMCMCbites( path = file.path(OutDir, outputs),
                                   burnin = 2,
                                   param.omit = c("s","z","sex","status"))
    
    
    ## Continue processing for this output inside the loop (e.g. extract mean values and CI for N, sigma etc...)
    # attach chains
    nimOutput <- do.call(rbind,nimOutput)
    # obtain columns names
    params <- colnames(nimOutput)
    # params.simple <- unique(sapply(strsplit(params, "\\["), "[", 1))
    # turn nimOutput as df 
    nimOutput <- as.data.frame(nimOutput)
    
    # run functions to get stats
    means<-col_mean(nimOutput, params)
    sd<-col_sd(nimOutput,params)
    CV<-col_cv(nimOutput,params)
    uci<-col_uci(nimOutput,params)
    lci<-col_lci(nimOutput,params)
    
    
    ## Store only the minimum in an object to summarize your results (e.g. N[sc,rp] <- mean(nimOutput[ ,"N"])
    # merge all stats in one df
    res <- as.data.frame(rbind(means,sd,CV,lci,uci))
    # creat columns for stats
    res2 <- tibble::rownames_to_column(res, "stat")
    params2 <- append(params, "stat", 0)
    colnames(res2) <- params2
    # store df per each repetition in a temporary list
    tempLista[[rp]] <- res2
    
  }#rp
  
  # store all scenarios in one final list
  resLista[[sc]] <- tempLista
  
}#sc


res25 <- do.call("rbind", resLista[[1]])
res25["scenario"] <- "25"
res50 <- do.call("rbind", resLista[[2]])
res50["scenario"] <- "50"
res75 <- do.call("rbind", resList[[3]])
res75["scenario"] <- "75"

all_res <- rbind(res25, res50,res75)
write.csv(all_res, "transect_sub.csv")


## -----------------------------------------------------------------------------
## ------   7. PLOT RESULTS -----

all_res <- read.csv(file.path(thisDir, "results_transectsub.csv"))
def_res <- read.csv(file.path(thisDir, "res5_2.csv"))
def_res["scenario"] <- "default"

all_res <- all_res[,-1]
def_res <- def_res[,-1]

def_colors <- ghibli::ghibli_palette("PonyoMedium")[4]

## ------   data wrangling

# means <- filter(all_res, stat == "means")
# sd <- filter(all_res, stat == "sd")
# CV <- filter(all_res, stat == "CV")
# lci <- filter(all_res, stat == "lci")
# uci <- filter(all_res, stat == "uci")


## -----------------------------------------------------------------------------
## ------   N   -----

# a <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5) +
#   geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) + 
#   labs(title="means", fill ="Repetitions simulation", y = "N", x = "Simulation") 
# 


g <- ggplot(data = all_res, aes(x=as.character(scenario), y=N, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  facet_wrap(vars(stat), scales = "free_y") + 
  geom_point(data = def_res, size =2, alpha = 1.2, color = def_colors) +
  # scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="N - Wolves abundance", fill ="Transects", y = "N", x = "Transects") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal") 

ggsave("transects_N_comp.jpeg", dpi = 300)

## -----------------------------------------------------------------------------
## ------   p0   -----

##------   data wrangling


p0 <- all_res[,13:18]
p0["stat"] <- all_res$stat
p0["scenario"] <- all_res$scenario
p0 <- melt(p0,id.vars=c("scenario","stat"))
# p01 <- p0 %>% 
#   filter(stat != "lci"| stat != "uci")


p0_52 <- def_res[,13:18]
p0_52["stat"] <- def_res$stat
p0_52["scenario"] <- def_res$scenario
p0_52 <- melt(p0_52,id.vars=c("scenario","stat"))


x_labs <- c("female alpha", "male alpha", "female pup",
            "male pup", "female other", "male other")


## ------   plots

# a <- ggplot(data = p0_means, aes(x=as.character(variable), y=value, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5)  +
#  scale_x_discrete(labels= x_labs) + labs(title="Means", x ="sex-status", y = "p0",fill = "Simulation")
# 


g <- ggplot(data = p0, aes(x=variable, y=value, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2, scale = "width") +
  geom_boxplot(width = 0.2, color="grey", alpha = 0.2) +
  geom_point(data = p0_52, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y") + 
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="p0 - baseline detectability", fill ="Search events", y = "p0", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_p0_comp.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   sigma   -----

## ------   data wrangling



sigma <- all_res[,21:26]
sigma["stat"] <- all_res$stat
sigma["scenario"] <- all_res$scenario
sigma <- melt(sigma,id.vars=c("scenario","stat"))

sigma_52 <- def_res[,21:26]
sigma_52["stat"] <- def_res$stat
sigma_52["scenario"] <- def_res$scenario
sigma_52 <- melt(sigma_52,id.vars=c("scenario","stat"))

## ------   plots

# a <- ggplot(data = sigma_means, aes(x=as.character(variable), y=value, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5) +
#   scale_x_discrete(labels= x_labs) +labs(title="Means", x ="sex-status", y = "sigma", fill = "Simulation")


g <- ggplot(data = sigma, aes(x=variable, y=value, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2, scale = "width") +
  geom_boxplot(width = 0.2, color="grey", alpha = 0.2) +
  geom_point(data = sigma_52, size =0.8, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y") + 
  scale_x_discrete(labels= x_labs) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="σ - scale parameter", fill ="Search events", y = "sigma", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_sigma_comp.jpeg", dpi = 300)

## -----------------------------------------------------------------------------
## ------   DetCov   -----

## ------   data wrangling

DetCov <- c("Transects Length", "Operator Experience", "Snow", "East/West", "Log Human Pop")

detcov <- all_res[,3:7]
detcov["stat"] <- all_res$stat
detcov["scenario"] <- all_res$scenario
detcov <- melt(detcov,id.vars=c("scenario","stat"))

detcov_52 <- def_res[,3:7]
detcov_52["stat"] <- def_res$stat
detcov_52["scenario"] <- def_res$scenario
detcov_52 <- melt(detcov_52,id.vars=c("scenario","stat"))

## ------   plots

# a <- ggplot(data = detcov_means, aes(x=as.character(variable), y=value, fill = as.character(scenario))) 
# a + geom_violin(alpha = 0.5) +
#   scale_x_discrete(labels= DetCov) +labs(title="Means", x ="sex-status", y = "Detection Covariates", fill = "Simulation")


g <- ggplot(data = detcov, aes(x=variable, y=value, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = detcov_52, size =0.8, alpha = 1.2, color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y") + 
  scale_x_discrete(labels= DetCov) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="Detection Covariates", fill ="Search events", y = "Detection Covariates", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_detcov_comp.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   HabCov   -----

## ------   data wrangling
HabCov <- c("Bare Rocks","Herbaceous", "Forest","Human Pop","Wolf presence")

habcov <- all_res[,8:12]
habcov["stat"] <- all_res$stat
habcov["scenario"] <- all_res$scenario
habcov <- melt(habcov,id.vars=c("scenario","stat"))

habcov_52 <- def_res[,8:12]
habcov_52["stat"] <- def_res$stat
habcov_52["scenario"] <- def_res$scenario
habcov_52 <- melt(habcov_52,id.vars=c("scenario","stat"))




g <- ggplot(data = habcov, aes(x=variable, y=value, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = habcov_52, size =0.8, alpha = 1.2, color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y") + 
  scale_x_discrete(labels= HabCov) +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="Detection Covariates", fill ="Search events", y = "Detection Covariates", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_habcov_comp.jpeg", dpi = 300)

## -----------------------------------------------------------------------------
## ------   psi   -----

g <- ggplot(data = all_res, aes(x=as.character(scenario), y=psi, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = def_res,size =2, alpha = 1.2,color = def_colors) +
  facet_wrap(vars(stat), scales = "free_y")  +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="ψ - real individual probability", fill ="Search events", y = "psi", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")


ggsave("transectsrepetition_psi_comp.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   rho   -----

g <- ggplot(data = all_res, aes(x=as.character(scenario), y=rho, fill = as.character(scenario))) 

g + geom_violin(alpha = 1.2) +
  geom_boxplot(width = 0.1, color="grey", alpha = 0.2) +
  geom_point(data = def_res, size =2, alpha = 1.2,color=def_colors) +
  facet_wrap(vars(stat), scales = "free_y") +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3)) +
  labs(title="ρ - sex proproportion", fill ="Search events", y = "rho", x = "Transects search events") +
  theme_ipsum() + 
  theme(plot.title = element_text(size=11),
        axis.text.x = element_text(size=8,hjust=1),
        axis.text.y = element_text(size=8),
        legend.position="bottom", legend.box = "horizontal")

ggsave("transectsrepetition_rho_comp.jpeg", dpi = 300)


## -----------------------------------------------------------------------------
## ------   8. FAILED ATTEMPTES-------------------------------------------------------------------

# load(file.path(thisDir, "parm.df.RData"))
# 
# inFiles <- list.files(path.IN)
# outFiles <- list.files(path.OUT)[grep("NimbleOutFOR", list.files(path.OUT))]
# ref <- expand.grid(list(G = c(F,T), H = c(F,T), A = c(F,T)))
# for(s in 1:nrow(parm.df)){
#   parm.df$scenario[s] <- which(ref$G == parm.df$Grouping[s] &
#                                  ref$H == parm.df$Heterogeneity[s] &
#                                  ref$A == parm.df$Adaptive[s])
# }#s


# for (p in colnames(nimOutput)) {

# means[p] <- apply(nimOutput[,p],2, mean)
# sd[p] <- apply(nimOutput[,p],2, sd)
# lci[p] <- apply(nimOutput[,p],2, quantile, probs = 0.05,  na.rm = TRUE)
# uci[p] <- apply(nimOutput[,p],2, quantile, probs = 0.95,  na.rm = TRUE)
# CV[p] <- apply(nimOutput[,p],2,FUN = mean()/sd())

# results[[p]][sc,rp] <- mean(nimOutput[, p])
# results[[p]][sc,rp] <- sd(nimOutput[, p])
# results[[p]][sc,rp] <- lci(nimOutput[, p])
# results[[p]][sc,rp] <- uci(nimOutput[, p])

N[st][rp,sc] <- mean(nimOutput[, "N"])
N[st][rp,sc] <- qs(nimOutput[, "N"],0.025)
N[st][rp,sc] <- qs(nimOutput[, "N"],0.975)
N[st][rp,sc] <- sd(nimOutput[, "N"])/mean(nimOutput[,"N"])
# p0RI_m[rp,sc] <- mean(nimOutput[, "p0[1, 2]"])
# p0RI_m[rp,sc] <- sd(nimOutput[, "p0[1, 2]"])
# p0RI_m[rp,sc] <- qs(nimOutput[, "p0[1, 2]"],0.025)
# p0RI_m[rp,sc] <- qs(nimOutput[, "p0[1, 2]"],0.975)
# p0RI_m[rp,sc] <- sd(nimOutput[, "p0[1, 2]"])/mean(nimOutput[,"N"])
# p0RI_f[rp,sc] <- mean(nimOutput[, "p0[1, 2]"])
# p0RI_f[rp,sc] <- sd(nimOutput[, "p0[1, 2]"])
# p0RI_f[rp,sc] <- qs(nimOutput[, "p0[1, 2]"],0.025)
# p0RI_f[rp,sc] <- qs(nimOutput[, "p0[1, 2]"],0.975)
# p0RI_f[rp,sc] <- sd(nimOutput[, "p0[1, 2]"])/mean(nimOutput[,"N"])

# }#p


# results <- list()
# for(k in 1:length(colnames.sims)){
#   results[[k]] <- length(rpp)
#   results[[k]]$mean <- rep_len(NA, length.out = rpp)
#   results[[k]]$lci <- rep_len(NA, length.out = rpp)
#   results[[k]]$uci <- rep_len(NA, length.out = rpp)
#   results[[k]]$sd <- rep_len(NA, length.out = rpp)
#   results[[k]]$CV <- rep_len(NA, length.out = rpp)
# }#k


# ## PROCESS OUTPUTS
myParams <- c("N","psi","rho", "p0", "theta", "sigma", "BetaHab", "BetaDet" )
# 
myResults <- list()
for(p in 1:length(myParams)){
  myResults[[p]] <- parm.df
  myResults[[p]]$mean <- rep(NA, nrow(parm.df))
  myResults[[p]]$lci <- rep(NA, nrow(parm.df))
  myResults[[p]]$uci <- rep(NA, nrow(parm.df))
  myResults[[p]]$sd <- rep(NA, nrow(parm.df))
  myResults[[p]]$CV <- rep(NA, nrow(parm.df))
  myResults[[p]]$Rhat <- rep(NA, nrow(parm.df))
}#p

n.rep <- max(parm.df$rep_ID)
n.set <- max(parm.df$set_ID)
for(s in 1:n.set){
  repCount <- 0
  for(r in 1:n.rep){
    fileName <- paste0("proc_NimbleOutFORNimbleInFile_Set",s,"_Rep",r,".RData")
    if(fileName %in% outFiles){
      load(file.path(path.OUT, fileName))
      if(all(na.omit(unlist(nimOutput$Rhat)) <= 1.1) & !is.na(any(unlist(nimOutput$Rhat) > 1.1)) ){
        repCount <- repCount + 1
        for(p in 1:length(myParams)){
          myResults[[p]]$mean[myResults[[p]]$set_ID == s & 
                                myResults[[p]]$rep_ID == r] <- nimOutput$mean[[myParams[p]]]
          myResults[[p]]$lci[myResults[[p]]$set_ID == s &
                               myResults[[p]]$rep_ID == r] <- nimOutput$q2.5[[myParams[p]]]
          myResults[[p]]$uci[myResults[[p]]$set_ID == s &
                               myResults[[p]]$rep_ID == r] <- nimOutput$q97.5[[myParams[p]]]
          myResults[[p]]$sd[myResults[[p]]$set_ID == s &
                              myResults[[p]]$rep_ID == r] <- nimOutput$sd[[myParams[p]]]
          myResults[[p]]$CV[myResults[[p]]$set_ID == s &
                              myResults[[p]]$rep_ID == r] <- nimOutput$sd[[myParams[p]]]/nimOutput$mean[[myParams[p]]]
          myResults[[p]]$Rhat[myResults[[p]]$set_ID == s &
                                myResults[[p]]$rep_ID == r] <- nimOutput$Rhat[[myParams[p]]]
          
        }#p
      }
    }#if
  }#r
  print(paste0("Set: ", s, " == ", repCount, " repetitions"))
}#s
save(myResults, file = file.path(thisDir, "results.RData"))



## ------   9. PIERRE'S HINTS -----



load(file.path(thisDir, "parm.df.RData"))
load(file.path(thisDir, "results.RData"))

## PLOTS N
{ 
  graphics.off()
  pdf(file = file.path(thisDir, "results_N.pdf"),
      width = 15, height = 8)
  
  species <- unique(parm.df$species)
  scenarios <- 1:8
  valNoise <- seq(-0.35,0.35, length.out = 8)
  par(mfrow = c(1,2))
  
  ## RELATIVE BIAS
  for(se in c("high", "low")){
    ylim <- 0.5
    plot(1,1, xlim = c(0.5,3.5), ylim = c(-ylim,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "RB(N)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(-ylim,ylim,length.out = 5),
         labels = seq(-ylim,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.2))
    abline(h = 0, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    
    p <- 1
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc  &
                                 myResults[[p]]$searchEffort == se, ]
        try(plot.violins2( dat.list = list(na.omit((temp$mean - temp$sim.N)/temp$sim.N)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 2.8, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  
  ## COEFFICIENT OF VARIATION
  for(se in c("high", "low")){
    ylim <- 0.15
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "CV(N)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.2))
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se,  ]
        try(plot.violins2( dat.list = list(na.omit(temp$CV)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 2.8, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  ## COVERAGE
  for(se in c("high", "low")){
    ylim <- 1
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "95%CI coverage", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.2))
    abline(h = 0.95, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se, ]
        
        coverage <- mean(temp$lci <= temp$sim.N & temp$uci >= temp$sim.N, na.rm =T)
        try(points( y = coverage,
                    x = sp + valNoise[sc],
                    pch = 21, cex = 2,
                    bg = colors[sc],
                    col = myCols[sc]), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 0.5, y = 0.35,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  graphics.off()
}

## PLOTS Sex Ratio
{ 
  graphics.off()
  pdf(file = file.path(thisDir, "results_SR.pdf"), width = 15, height = 8)
  
  species <- unique(parm.df$species)
  scenarios <- 1:8
  valNoise <- seq(-0.35,0.35, length.out = 8)
  par(mfrow = c(1,2))
  
  ## RELATIVE BIAS
  for(se in c("high", "low")){
    ylim <- 0.5
    plot(1,1, xlim = c(0.5,3.5), ylim = c(-ylim,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "RB(Sex Ratio)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(-ylim,ylim,length.out = 5),
         labels = seq(-ylim,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.1))
    abline(h = 0, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    p <- 3
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc  &
                                 myResults[[p]]$searchEffort == se, ]
        try(plot.violins2( dat.list = list(na.omit((temp$mean - temp$sim.SR)/temp$sim.SR)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 2.8, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  
  ## COEFFICIENT OF VARIATION
  for(se in c("high", "low")){
    ylim <- 0.15
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "CV(N)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.1))
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se,  ]
        try(plot.violins2( dat.list = list(na.omit(temp$CV)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 3, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  ## COVERAGE
  for(se in c("high", "low")){
    ylim <- 1
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "95%CI coverage", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.1))
    abline(h = 0.95, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se, ]
        
        coverage <- mean(temp$lci <= temp$sim.SR & temp$uci >= temp$sim.SR, na.rm = T)
        try(points( y = coverage,
                    x = sp + valNoise[sc],
                    pch = 21, cex = 2,
                    bg = colors[sc],
                    col = myCols[sc]), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 0.5, y = 0.35,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  graphics.off()
}



## -----------------------------------------------------------------------------
