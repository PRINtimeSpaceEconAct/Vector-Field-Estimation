rm(list = ls())
source("loadLib.R")


nObs = 100000
X = matrix(nrow=nObs,rnorm(2*nObs))
nEval = 2500
xGrid = seq(from=min(X[,1]),to=max(X[,1]),length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X[,2]),to=max(X[,2]),length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid,yGrid))


system.time(est <- densityEst2d(X,x,kernel = "gauss", h = 0.05, gc=TRUE))

# Create a 3D scatter plot
xCoord = est$x[,1]
yCoord = est$x[,2]
z = est$densityEst
plot_ly(x = xCoord,y = yCoord,z = z, intensity = z,type = "mesh3d") 

library(sm)
system.time(est.sm <- sm.density(X,h=c(0.05,0.05),eval.points = x,eval.grid=FALSE,nbins=0))
z = est.sm$estimate-est$densityEst
plot_ly(x = xCoord,y = yCoord,z = z, intensity = z,type = "mesh3d") 

                    