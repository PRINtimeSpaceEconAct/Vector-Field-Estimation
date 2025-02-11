# Clear workspace and load dependencies
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
rm(list = ls())
DEBUG = TRUE
source("src/libs/loadLib.R")
library(dplyr)
library(sm)

load("datasets/datasetGDP.LE.NonOverlapping.RData")
# load("datasets/datasetGDP.LE.Overlapping.RData")
data = dataset.GDP.LE
X0 = cbind(as.numeric(data$GDP.t0),as.numeric(data$LE.t0))
X1 = cbind(as.numeric(data$GDP.t1),as.numeric(data$LE.t1))

# data = dataset.GDP.LE.REL
# X0 = cbind(as.numeric(data$GDP.REL.t0),as.numeric(data$LE.REL.t0))
# X1 = cbind(as.numeric(data$GDP.REL.t1),as.numeric(data$LE.REL.t1))

# draw transitions as they are ----
plot(X0, type = "n", xlab = "X", ylab = "Y", main = "True Vector Field")
arrows(X0[,1], X0[,2], X1[,1], X1[,2], length = 0.05, angle = 15, col = "black")

# parameters ----
nEval = 2500

# eval points
xGrid = seq(from=-20000, to=1.1*max(c(X0[,1],X1[,1])), length.out=round(sqrt(nEval)))
yGrid = seq(from=0.9*min(c(X0[,2],X1[,2])), to=1.1*max(c(X0[,2],X1[,2])), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))


# stima
est_field = LLfield(X0, X1, x=x, kernel.type="epa",method.h = "sj",
                                     chunk_size=1000,
                                     sparse=FALSE, gc=TRUE)


# plots ----
## plot campo stimato ----
lengthArrows = 1.0
latestYear = max(data$Year.end)
firstYear = min(data$Year.ini)
dataLatest = filter(data,Year.end == latestYear)
dataFirst = filter(data,Year.ini == firstYear)


dev.new()
# campo
plot(est_field$x, type = "n", xlab = "GDP per capita", ylab = "Life Expectancy", main = "Estimated Vector Field")
arrows(est_field$x[,1], est_field$x[,2],
       est_field$x[,1] + lengthArrows*est_field$estimator[,1], 
       est_field$x[,2] + lengthArrows*est_field$estimator[,2],
       length = 0.05, angle = 15, col = "blue")
# punti
points(dataFirst$GDP.t0,dataFirst$LE.t0,pch=19,col="black",cex=0.5)
# text(dataFirst$GDP.t0,dataFirst$LE.t0,labels=dataFirst$countryCode,cex=1,pos=4,col="black")
points(dataLatest$GDP.t0,dataLatest$LE.t0,pch=19,col="red",cex=0.5)
text(dataLatest$GDP.t0,dataLatest$LE.t0,labels=dataLatest$countryCode,cex=1,pos=4,col="red")

# non parametriche
estFirst = sm.regression(dataFirst$GDP.t0,dataFirst$LE.t0,display="none")
estLatest = sm.regression(dataLatest$GDP.t1,dataLatest$LE.t1,display="none")
lines(estFirst$eval.points,estFirst$estimate,col="black",lwd=4)
lines(estLatest$eval.points,estLatest$estimate,col="red",lwd=4)


# forecast delle obs 
speedFactor = 0.1
nPeriods = (2015-1960)/5 * 1/speedFactor
forecastObs = forecastDiscrete(cbind(dataFirst$GDP.t0,dataFirst$LE.t0), est_field, speedFactor = speedFactor, nPeriods = nPeriods)
points(forecastObs[,,nPeriods],col="purple",pch=19,cex=0.5)
estForecastObs = sm.regression(forecastObs[,1,nPeriods],forecastObs[,2,nPeriods],display="none")
lines(estForecastObs$eval.points,estForecastObs$estimate,col="purple",lwd=4)

# forecast della est 
forecastEst = forecastDiscrete(cbind(estFirst$eval.points,estFirst$estimate), est_field, speedFactor = speedFactor, nPeriods = nPeriods)
lines(forecastEst[,1,nPeriods],forecastEst[,2,nPeriods],col="green",lwd=4)

legend("bottomright",legend=c("1960","2015","Forecasted Observed","Forecasted Estimated"),
       col=c("black","red","purple","green"),pch=c(19,19,19,19),lwd=c(4,4,4,4),bg="white")
dev.copy2pdf(file="Forecast19602015NonOverlapping.pdf")

## plot direzioni verticali orizzontali campo vettoriale ----
library(fields)

Fx = est_field$estimator[,1]
dev.new()
image.plot(x = unique(est_field$x[,1]), y = unique(est_field$x[,2]), z = matrix(Fx, nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="GDP per capita",ylab="Life Expectancy",main="Fx")
points(dataFirst$GDP.t0,dataFirst$LE.t0,pch=19,col="black",cex=1.0)
points(dataLatest$GDP.t0,dataLatest$LE.t0,pch=23,col="black",cex=1.0)
text(dataLatest$GDP.t0,dataLatest$LE.t0,labels=dataLatest$countryCode,cex=1,pos=4,col="black")
estFirst = sm.regression(dataFirst$GDP.t0,dataFirst$LE.t0,display="none")
lines(estFirst$eval.points,estFirst$estimate,col="black",lwd=4)
estLatest = sm.regression(dataLatest$GDP.t1,dataLatest$LE.t1,display="none")
lines(estLatest$eval.points,estLatest$estimate,col="black",lwd=4,lty = 3)
dev.copy2pdf(file="Fx19602015NonOverlapping.pdf")

Fy = est_field$estimator[,2]
dev.new()
image.plot(x = unique(est_field$x[,1]), y = unique(est_field$x[,2]), z = matrix(Fy, nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="GDP per capita",ylab="Life Expectancy",main="Fy")
points(dataFirst$GDP.t0,dataFirst$LE.t0,pch=19,col="black",cex=1.0)
points(dataLatest$GDP.t0,dataLatest$LE.t0,pch=23,col="black",cex=1.0)
text(dataLatest$GDP.t0,dataLatest$LE.t0,labels=dataLatest$countryCode,cex=1,pos=4,col="black")
estFirst = sm.regression(dataFirst$GDP.t0,dataFirst$LE.t0,display="none")
lines(estFirst$eval.points,estFirst$estimate,col="black",lwd=4)
estLatest = sm.regression(dataLatest$GDP.t1,dataLatest$LE.t1,display="none")
lines(estLatest$eval.points,estLatest$estimate,col="black",lwd=4,lty = 3)
dev.copy2pdf(file="Fy19602015NonOverlapping.pdf")

FyOverFx = est_field$estimator[,2]/est_field$estimator[,1]
dev.new()
image.plot(x = unique(est_field$x[,1]), y = unique(est_field$x[,2]), z = matrix(FyOverFx, nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="GDP per capita",ylab="Life Expectancy",main="Fy/Fx",breaks=seq(from=quantile(FyOverFx,0.05),to=quantile(FyOverFx,0.95),length.out=65))
points(dataFirst$GDP.t0,dataFirst$LE.t0,pch=19,col="black",cex=1.0)
points(dataLatest$GDP.t0,dataLatest$LE.t0,pch=23,col="black",cex=1.0)
text(dataLatest$GDP.t0,dataLatest$LE.t0,labels=dataLatest$countryCode,cex=1,pos=4,col="black")
estFirst = sm.regression(dataFirst$GDP.t0,dataFirst$LE.t0,display="none")
lines(estFirst$eval.points,estFirst$estimate,col="black",lwd=4)
estLatest = sm.regression(dataLatest$GDP.t1,dataLatest$LE.t1,display="none")
lines(estLatest$eval.points,estLatest$estimate,col="black",lwd=4,lty = 3)
dev.copy2pdf(file="FyoverFx19602015NonOverlapping.pdf")

# VAR, questo è fatto a culo ----
library(vars)
dataEnd = filter(data,Year.end == max(data$Year.end))  %>% select(GDP = GDP.t1,LE = LE.t1, t = Year.end, countryCode)
dataVAR = select(data, GDP = GDP.t0, LE = LE.t0, t = Year.ini, countryCode)
dataVAR = rbind(dataVAR,dataEnd)

var_model <- VAR(select(dataVAR,c(GDP,LE)), p = 1, type = "const")  # p = 1 for VAR(1)
summary(var_model)

# Extract coefficients
coeffs <- coef(var_model)
C <- sapply(coeffs, function(x) x["const", 1])
A1 <- do.call(rbind, lapply(coeffs, function(mat) mat[grep("l1", rownames(mat)), 1, drop = TRUE]))

# VAR a mano, ma forse giusto ----
lm1 = lm(GDP.t1 ~ GDP.t0 + LE.t0, data = data)
lm2 = lm(LE.t1 ~ GDP.t0 + LE.t0, data = data)
A = matrix(c(coef(lm1)[2:3],coef(lm2)[2:3]),nrow=2,byrow=TRUE)
C = c(coef(lm1)[1],coef(lm2)[1])

# lm1 = lm(GDP.t1 ~ GDP.t0 + LE.t0 + 0, data = data)
# lm2 = lm(LE.t1 ~ GDP.t0 + LE.t0 + 0, data = data)
# A = matrix(c(coef(lm1)[1:2],coef(lm2)[1:2]),nrow=2,byrow=TRUE)
# C = c(0,0)

VFVAR = function(x){
  return(A %*% x + C)
}
VFVARx = t(apply(x, 1, VFVAR))



## plot ----
latestYear = max(data$Year.end)
firstYear = min(data$Year.ini)
dataLatest = filter(data,Year.end == latestYear)
dataFirst = filter(data,Year.ini == firstYear)

dev.new()
# campo
lengthArrowsVAR = 0.1
plot(est_field$x, type = "n", xlab = "GDP per capita", ylab = "Life Expectancy", main = "Estimated VAR deg 1 (length 0.01)")
arrows(est_field$x[,1], est_field$x[,2],
       est_field$x[,1] +  lengthArrowsVAR*VFVARx[,1], 
       est_field$x[,2] +  lengthArrowsVAR*VFVARx[,2],
       length = 0.05, angle = 15, col = "blue")

# punti
points(dataFirst$GDP.t0,dataFirst$LE.t0,pch=19,col="black",cex=0.5)
# text(dataFirst$GDP.t0,dataFirst$LE.t0,labels=dataFirst$countryCode,cex=1,pos=4,col="black")
points(dataLatest$GDP.t0,dataLatest$LE.t0,pch=19,col="red",cex=0.5)
text(dataLatest$GDP.t0,dataLatest$LE.t0,labels=dataLatest$countryCode,cex=1,pos=4,col="red")

# non parametriche
estFirst = sm.regression(dataFirst$GDP.t0,dataFirst$LE.t0,display="none")
estLatest = sm.regression(dataLatest$GDP.t1,dataLatest$LE.t1,display="none")
lines(estFirst$eval.points,estFirst$estimate,col="black",lwd=4)
lines(estLatest$eval.points,estLatest$estimate,col="red",lwd=4)


# forecast delle obs 
speedFactor = 0.1
nPeriods = (2015-1960)/5 * 1/speedFactor * 0.05
forecastObs = cbind(dataFirst$GDP.t0,dataFirst$LE.t0)
for (i in 1:nPeriods){
    forecastObs = forecastObs + speedFactor * t(apply(forecastObs, 1, VFVAR))
}

points(forecastObs,col="purple",pch=19,cex=0.5)
estForecastObs = sm.regression(forecastObs[,1],forecastObs[,2],display="none")
lines(estForecastObs$eval.points,estForecastObs$estimate,col="purple",lwd=4)

# forecast della est 
forecastEst = cbind(estFirst$eval.points,estFirst$estimate)
for (i in 1:nPeriods){
    forecastEst = forecastEst + speedFactor * t(apply(forecastEst, 1, VFVAR))
}
lines(forecastEst[,1],forecastEst[,2],col="green",lwd=4)

legend("bottomright",legend=c("1960","2015","Forecasted Observed for 1/20th","Forecasted Estimated for 1/20th"),
       col=c("black","red","purple","green"),pch=c(19,19,19,19),lwd=c(4,4,4,4),bg="white")
dev.copy2pdf(file="ForecastVAR19602015NonOverlapping.pdf")


# VAR con quadrato a mano, ma forse giusto ----
lm1 = lm(GDP.t1 ~ GDP.t0 + LE.t0 + I(GDP.t0^2) + I(LE.t0^2), data = data)
lm2 = lm(LE.t1 ~ GDP.t0 + LE.t0 + I(GDP.t0^2) + I(LE.t0^2), data = data)

A1 = matrix(c(coef(lm1)[2:3],coef(lm2)[2:3]),nrow=2,byrow=TRUE)
A2 = matrix(c(coef(lm1)[4:5],coef(lm2)[4:5]),nrow=2,byrow=TRUE)
C = c(coef(lm1)[1],coef(lm2)[1])

VFVAR = function(x){
    return(A1 %*% x + A2 %*% x^2 + C)
}
VFVARx = t(apply(x, 1, VFVAR))


lengthArrowsVAR = 0.1
plot(est_field$x, type = "n", xlab = "GDP per capita", ylab = "Life Expectancy", main = "Estimated VAR deg 2")
arrows(est_field$x[,1], est_field$x[,2],
       est_field$x[,1] + lengthArrowsVAR*VFVARx[,1], 
       est_field$x[,2] + lengthArrowsVAR*VFVARx[,2],
       length = 0.05, angle = 15, col = "blue")




