#==============================================================================#
#                                                                              #
#                 Occupancy data improves parameter precision in               #
#                       spatial capture-recapture models                       #
#                    ~ Stone marten in Cabañeros National Park ~               #
#                              ECOLOGY AND EVOLUTION                           #
#                José Jiménez, Francisco Díaz-Ruiz, Pedro Monterroso,          #
#                        Jorge Tobajas, Pablo Ferreras                         #
#                             9/08/2022 10:33:56 a. m.                         #
#                                                                              #
#==============================================================================#

setwd('...') # Define your work directory, with uncompressed data, and 
             # auxiliary functions (Funciones__SCR.R and nimSummary.R)
source("Funciones_SCR.R")
source("nimSummary.R")
library(lattice)
library(secr)


# CAMERA TRAPS
#==================
Oper.cam<-data.matrix(read.table("Oper_cam.txt", header=FALSE)); dim(Oper.cam)
KT.cam<-apply(Oper.cam,1,sum)

cam<-data.matrix(read.table("cam_traps.txt", header=FALSE))
cam<-cam[,2:3]
head(cam)
X.cam<-  cam/1000
X1.mean<- mean(X.cam[,1])
X2.mean<- mean(X.cam[,2])

X.cam[,1]<-X.cam[,1]-X1.mean
X.cam[,2]<-X.cam[,2]-X2.mean

J.cam<-nrow(X.cam)
K.cam<-ncol(Oper.cam)

# Camera trap operation
J.trap<-60
x1<-as.matrix(Oper.cam)  ## Tenemos que convertir esto en matriz
image(1:K.cam,1:J.cam,t(x1), yaxt = "n", xlab="Occasion", ylab="", cex.lab=1.25, col=topo.colors(2))
mtext(side = 2, "Camera-trap", line = 2.5, cex=1.25)
axis(2, rev(seq(1, J.trap, by=2)))
box()


# BOX-TRAPS
#==================
# Trapping occasions: 121 days (15/01/2014-22/04/2014)
Oper.trap<-data.matrix(read.table("Oper_trap.txt", header=FALSE))
dim(Oper.trap)

KT.trap<-apply(Oper.trap,1,sum)


# BOX TRAP DATA
# We used secr (Efford, 2022) to prepare and format data:
StoneMarten.ch <- read.capthist("capt.txt", "box_traps.txt", detector='count', noccasions=121)
summary(StoneMarten.ch)

traplocs<-as.matrix(traps(StoneMarten.ch))
X.trap<-data.matrix(traplocs)/1000
X.trap[,1]<-(X.trap[,1]-X1.mean)
X.trap[,2]<-(X.trap[,2]-X2.mean)

StoneMarten<-aperm(StoneMarten.ch, c(1,3,2))
str(StoneMarten)

datYknown<- as.array(StoneMarten)

cat("Box trap - days=", sum(Oper.trap), "\n") # box trap-day
cat("Mean trapping days (box trap)=", mean(apply((Oper.trap),1,sum)), "\n")

K.trap<-dim(Oper.trap)[2] ;K.trap     # box trap sampling occasions
J.trap<-dim(Oper.trap)[1] ;J.trap     # number of camera traps

# Box trap operation
x<-as.matrix(Oper.trap)  ## in matrix format
image(1:K.trap,1:J.trap,t(x), yaxt = "n", xlab="Occasion", ylab="", cex.lab=1.25, col=topo.colors(2))
mtext(side = 2, "Live trap", line = 2.5, cex=1.25)
axis(2, rev(seq(1, J.trap, by=2)))
box()

# Operation summarised by sampling occasion
KT.trap<-apply(Oper.trap,1,sum)
colnames(Oper.trap)<-c(1:K.trap)

M<-500 # data augmentation
yr<-apply(datYknown, c(2,3),sum) # To test captures by operation
which(yr*Oper.trap -yr<0, arr.in=TRUE)

nind<- dim(datYknown)[1]; nind     # captured and marked individuals
Yaug <- array(0, dim = c(M, J.trap, K.trap))
Yaug[1:nind, , ] <- datYknown
y<-apply(Yaug, 1:2, sum)


# Capture plot
plot(X.cam, pch="+", col="blue", xlim=c(min(X.cam[,1])-1,max(X.cam[,1])+1), 
  ylim=c(min(X.cam[,2])-1,max(X.cam[,2])+1),  
  xlab="X",
  ylab="Y",
  asp=TRUE)
points(X.trap, pch="+", col="red")
# spiderplot
n1<-apply(datYknown,2,sum)
symbols(X.trap, circles=n1/5, inches=F,bg="#0000FF3F", fg=NULL, add=T)


# CAMERA TRAP DATA
## UNMARKED
(K.cam<-dim(Oper.cam)[2])
n3d<-secr::read.capthist("detections.txt", "cam_traps.txt", noccasions =99)
nb<-apply(n3d,c(3,2),sum)
# We trasnformed counts to binary detections:
nb[nb>1]<-1
# To test captures by operation
which(nb*Oper.cam -nb<0, arr.in=TRUE)
nred<-apply(nb,1,sum)

# Detection plot
tot2<-apply(nb,1,sum)
symbols(X.cam, circles=tot2/5, inches=F,bg="#EEAD0E66", fg=NULL, add=T)
# Spiderplot
spiderplotJJ4(datYknown, X.trap, lwd=2,buffer=1)
 
# TELEMETRY
#===============
library(adehabitatHR)
locs<-read.table("PosStoneMarten.txt", header=TRUE)[,2:3]/1000
locs[,1]<-locs[,1]-X1.mean
locs[,2]<-locs[,2]-X2.mean
locs1<-locs[1:12,]
locs2<-locs[13:21,]
locs3<-locs[22:42,]
locs4<-locs[43:56,]

points(locs1[,1], locs1[,2], pch=16, col="#FF00004C", cex=1)
points(locs2[,1], locs2[,2], pch=16, col="#0000FF4C", cex=1)

xy1<-cbind(locs1[,1], locs1[,2])
xy1<-data.matrix(xy1)
xysp1 <- SpatialPoints(xy1)
kudl1 <- kernelUD(xysp1, grid=100)
homerange1 <- getverticeshr(kudl1); homerange1@data$area*1E6

xy2<-cbind(locs2[,1], locs2[,2])
xy2<-data.matrix(xy2)
xysp2 <- SpatialPoints(xy2)
kudl2 <- kernelUD(xysp2, grid=100)
homerange2 <- getverticeshr(kudl2); homerange2@data$area*1E6
#contour(getvolumeUD(kudl2)[1], level=c(50,95), add=TRUE, col=c('blue','blue'), lwd=2)

xy3<-cbind(locs3[,1], locs3[,2])
xy3<-data.matrix(xy3)
xysp3 <- SpatialPoints(xy3)
kudl3 <- kernelUD(xysp3, grid=100)
homerange3 <- getverticeshr(kudl3); homerange3@data$area*1E6
#contour(getvolumeUD(kudl2)[1], level=c(50,95), add=TRUE, col=c('blue','blue'), lwd=2)

xy4<-cbind(locs4[,1], locs4[,2])
xy4<-data.matrix(xy4)
xysp4 <- SpatialPoints(xy4)
kudl4 <- kernelUD(xysp4, grid=100)
homerange4 <- getverticeshr(kudl4); homerange4@data$area*1E6
#contour(getvolumeUD(kudl2)[1], level=c(50,95), add=TRUE, col=c('blue','blue'), lwd=2)

points(locs3[,1], locs3[,2], pch=16, col="#00FF004C", cex=1)
points(locs4[,1], locs4[,2], pch=16, col="forestgreen", cex=1)

locs<-rbind(locs1,locs2,locs3,locs4)
inds<-c(rep(1,nrow(locs1)),rep(2,nrow(locs2)),rep(3,nrow(locs3)),rep(4,nrow(locs4)))
nlocs<-nrow(locs)
n.collar<-4

buff<- 2 # ( buff>3*sigma )
xl<-min(X.cam[,1])-buff
xu<-max(X.cam[,1])+buff
yl<-min(X.cam[,2])-buff
yu<-max(X.cam[,2])+buff
xlims=c(xl, xu)
ylims=c(yl, yu)
area<-(xu-xl)*(yu-yl)



library(nimble)
## define the model using NIMBLE (de Valpine et al., 2017)
code <- nimbleCode({

  p0.trp ~ dunif(0,1)
  p0.cam ~ dunif(0,5)
  sigma ~ dunif(0,5)
  psi ~ dunif(0,1)

  # Capture data and euclidean distances
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    dtrap[i,1:J.trap] <- (s[i,1]-X.trap[1:J.trap,1])^2 + (s[i,2]-X.trap[1:J.trap,2])^2
    p[i,1:J.trap] <-p0.trp* exp(- dtrap[i,1:J.trap]/(2*sigma^2))
    dcam[i,1:J.cam] <- (s[i,1] - X.cam[1:J.cam,1])^2 + (s[i,2] - X.cam[1:J.cam,2])^2
    lam[i,1:J.cam] <- p0.cam*exp(-dcam[i,1:J.cam]/(2*sigma^2))

    # Sub-model for capture-recapture
    for(j in 1:J.trap){
      pscr[i,j]  <- 1 - exp(-p[i,j])
      y[i,j] ~ dbinom(pscr[i,j]*z[i], KT.trap[j])
    }
    # Compute detection probability for occupancy
    for(j in 1:J.cam) {
      pocc[i,j] <- 1 - exp(-lam[i,j])
      # for PA data compute probability of not captured
      pn[i,j] <- (1 - (pocc[i,j]*z[i]) )
    } #j
  }
  # Sub-model for the presence-absence data
  for(j in 1:J.cam) {
    yocc[j] ~ dbinom((1-prod(pn[1:M,j])),KT.cam[j]) #
  }#j 

  # Sub-model for telemetry data
  for (r in 1:nlocs){
    locs[r,1]~dnorm(s[inds[r],1], 1/(sigma^2))
    locs[r,2]~dnorm(s[inds[r],2], 1/(sigma^2))
  }

  N <- sum(z[1:M])
  D <- N/area
})

str(constants <- list(M = M,
                      J.trap=J.trap,
                      J.cam=J.cam,
                      area=area,
                      K.trap=121,
                      K.cam=99,
                      nlocs=nlocs,
                      inds=inds))

str ( data  <-  list (y=y,
                      yocc=nred,
                      X.trap=X.trap,
                      X.cam=X.cam,
                      KT.trap=KT.trap,
                      KT.cam=KT.cam,
                      xlim=xlims,
                      ylim=ylims,
                      locs=locs))

## Init for s
sst <- cbind(runif(M,xlims[1],xlims[2]),runif(M,ylims[1],ylims[2]))
for(i in 1:nind){
  sst[i,1] <- mean( X.trap[y[i,]>0,1] )
  sst[i,2] <- mean( X.trap[y[i,]>0,2] )
}
zst<-rep(1,M)

str ( inits   <- list (p0.trp=0.1,
                       p0.cam=0.05,
                       sigma=2,
                       psi=0.3,
                       s=sst,
                       z=zst))

params<- c('p0.trp','p0.cam','psi','N','D','sigma')

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
Rmodel$initializeInfo()
#Rmodel$calculate()
Cmodel <- compileNimble(Rmodel)

conf<-configureMCMC(Rmodel, monitors = params, useConjugacy=TRUE, thin=1)

conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'AF_slice',control=list(adaptive=TRUE,adaptScaleOnly=TRUE),silent = TRUE)
}

conf$removeSamplers('z')
for(node in Rmodel$expandNodeNames('z')) conf$addSampler(target = node, type = 'slice')

MCMC <- buildMCMC(conf)
Cmcmc <- compileNimble(MCMC, project = Rmodel)

nb = 10000
ni = 50000 + nb
nc = 3

outNim <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc,
                  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

nimSummary(outNim,digits=3) # Rhat is the potential scale reduction factor on 
                            # rank normalized split chains (Vehtari et al., 2020)


# REFERENCES
#=============
# de Valpine, P., D. Turek, C.J. Paciorek, C. Anderson-Bergman, D.
#   Temple Lang, and R. Bodik. 2017. Programming with models: writing
#   statistical algorithms for general model structures with NIMBLE.
#   Journal of Computational and Graphical Statistics 26: 403-413.
#   <DOI:10.1080/10618600.2016.1172487>.
  
# Efford, M. G. (2022). secr: Spatially explicit capture-recapture
#   models. R package version 4.5.5.
#   https://CRAN.R-project.org/package=secr

# Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2020). 
#   Rank-Normalization, Folding, and Localization: An Improved R for Assessing 
#   Convergence of MCMC (with Discussion). Bayesian Analysis, 16(2), 667-718. 
#   https://doi.org/10.1214/20-ba1221