# Clear workspace and load dependencies
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
# setwd("/Users/davidefiaschi/Library/CloudStorage/OneDrive-UniversityofPisa/Cristiano\ Ricci\'s\ files\ -\ timeSpaceEvolutionEcAct/RVF/R\ code/Vector\ Field\ Estimation")

rm(list = ls())
DEBUG = FALSE
source("src/libs/loadLib.R")
library(dplyr)
library(sm)
library(xtable)

timeInterval = 5 # 5 or 10 years
load("datasets/datasetGDP.LE.NonOverlapping.RData")
# load("datasets/datasetGDP.LE.Overlapping.RData")
data = dataset.GDP.LE

if (timeInterval == 10) {
  df1 <- data
  df2 <- data
  
  data_10y <- inner_join(
    df1, 
    df2, 
    by = "countryCode",
    suffix = c(".x", ".y"),
    relationship = "many-to-many"
  ) %>%
    filter(Year.end.x == Year.ini.y) %>%
    transmute(
      countryCode = countryCode,
      Year.ini = Year.ini.x,
      Year.end = Year.end.y,
      GDP.t0 = GDP.t0.x,
      LE.t0 = LE.t0.x,
      GDP.t1 = GDP.t1.y,
      LE.t1 = LE.t1.y
    )
  
  # Filter to get non-overlapping 10-year intervals (e.g., 1960-1970, 1970-1980)
  startYearOfDecade = min(data_10y$Year.ini) %% 10
  data = data_10y %>% filter(Year.ini %% 10 == startYearOfDecade)
}

countriesDavide = read.csv(file = "src/countryCodeAnalysis.csv")

# Filter countries in Davide's file
data = data %>% filter(countryCode %in% countriesDavide[,2])

# count countries
length(unique(data$countryCode))

# data = data %>% mutate(GDP.t0 = log(GDP.t0), GDP.t1 = log(GDP.t1))
summary(log(data$GDP.t1/data$GDP.t0))
X0 = cbind(as.numeric(data$GDP.t0),as.numeric(data$LE.t0))
X1 = cbind(as.numeric(data$GDP.t1),as.numeric(data$LE.t1))

# tabella country code
# codiciAlpha3 = read.csv("datasets/codici ALPHA3.csv")
# data = data %>% left_join(codiciAlpha3, by = c("countryCode" = "Code"))
# # remove leading blank space from country names
# data$Country = gsub("^\\s+|\\s+$", "", data$Country)
# MCodes = as.matrix(unique(select(data, c(countryCode, Country))))
# 
# MCodesLarge = matrix("",nrow=40,ncol = 6)
# MCodesLarge[1:40,1:2] = MCodes[1:40,]
# MCodesLarge[1:40,3:4] = MCodes[41:80,]
# MCodesLarge[1:25,5:6] = MCodes[81:105,]
# colnames(MCodesLarge) = c("Code","Country","Code","Country","Code","Country")
# xtable(MCodesLarge)


# data = dataset.GDP.LE.REL
# X0 = cbind(as.numeric(data$GDP.REL.t0),as.numeric(data$LE.REL.t0))
# X1 = cbind(as.numeric(data$GDP.REL.t1),as.numeric(data$LE.REL.t1))

# draw transitions as they are ----
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
dataFirst = filter(data,Year.ini == min(data$Year.ini))
plot(dataFirst$GDP.t0,dataFirst$LE.t0, type = "n", xlab = "GDP per capita (PPP in million 2017 USD)", ylab = "Life expectancy at birth", main = "",xlim=range(dataFirst$GDP.t0,dataFirst$GDP.t1),ylim=range(dataFirst$LE.t0,dataFirst$LE.t1))
arrows(dataFirst$GDP.t0,dataFirst$LE.t0, dataFirst$GDP.t1,dataFirst$LE.t1, length = 0.15, angle = 15, col = "black")
points(dataFirst$GDP.t0,dataFirst$LE.t0,pch=19,col="black",cex=0.5)
points(dataFirst$GDP.t1,dataFirst$LE.t1,pch=19,col="red",cex=0.5)
grid()
#dev.copy2pdf(file="scatterArrowsPreston19601965.pdf")

# parameters ----
nEval = 2500

# eval points
# xGrid = seq(from=-20000, to=1.1*max(c(X0[,1],X1[,1])), length.out=round(sqrt(nEval)))
# yGrid = seq(from=0.9*min(c(X0[,2],X1[,2])), to=1.2*max(c(X0[,2],X1[,2])), length.out=round(sqrt(nEval)))
# x = as.matrix(expand.grid(xGrid, yGrid))

xGrid = seq(from=0.9*min(c(X0[,1],X1[,1])), to=1.1*max(c(X0[,1],X1[,1])), length.out=round(sqrt(nEval)))
yGrid = seq(from=0.9*min(c(X0[,2],X1[,2])), to=1.1*max(c(X0[,2],X1[,2])), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))


# stima
# est_field = LLfield(X0, X1, x=x, kernel.type="epa",method.h = "sj",
#                                      chunk_size=1000,
#                                       gc=TRUE)
# est_field = NWfield(X0, X1, x=x, kernel.type="epa",h=1.0,
#                                         chunk_size=1000,
#                                          gc=TRUE
#                     #,method.h="sj"
#                     )
est_field = NWfield(X0, X1, x=x, kernel.type="epa",h=0.6418,
                                        chunk_size=1000,
                                        gc=TRUE)


# plots ----
## plot campo stimato ----
lengthArrows = 1.0
latestYear = max(data$Year.end)
firstYear = min(data$Year.ini)
dataLatest = filter(data,Year.end == latestYear)
dataFirst = filter(data,Year.ini == firstYear)


dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 

# campo
#signifVFest = significanceVF(est_field,X0,X1,0.05)
signifVFest = significanceVF(est_field,p_crit=0.05)

plot(est_field$x, type = "n", xlab = "GDP per capita (PPP in million 2017 USD, log scale)", ylab = "Life Expectancy", main = "",xaxt="n")
axis(1, at=6:12, labels=as.character(round(exp(6:12),-3)))
arrows(est_field$x[,1], est_field$x[,2],
       est_field$x[,1] + lengthArrows*signifVFest$signif*est_field$estimator[,1], 
       est_field$x[,2] + lengthArrows*signifVFest$signif*est_field$estimator[,2],
       length = 0.05, angle = 15, col = "black",lwd=1)
# punti
points(dataFirst$GDP.t0,dataFirst$LE.t0,pch=19,col="blue",cex=0.5)
# text(dataFirst$GDP.t0,dataFirst$LE.t0,labels=dataFirst$countryCode,cex=1,pos=4,col="black")
points(dataLatest$GDP.t0,dataLatest$LE.t0,pch=19,col="red",cex=0.5)
text(dataLatest$GDP.t0,dataLatest$LE.t0,labels=dataLatest$countryCode,cex=0.5,pos=4,col="red")

# non parametriche
estFirst = sm.regression(dataFirst$GDP.t0,dataFirst$LE.t0,display="none")
estLatest = sm.regression(dataLatest$GDP.t1,dataLatest$LE.t1,display="none")
lines(estFirst$eval.points,estFirst$estimate,col="blue",lwd=3)
lines(estLatest$eval.points,estLatest$estimate,col="red",lwd=3)


# forecast delle obs 
speedFactor = 0.1
nPeriods = (2015-1960)/5 * 1/speedFactor
forecastObs = forecastDiscrete(cbind(dataFirst$GDP.t0,dataFirst$LE.t0), est_field, speedFactor = speedFactor, nPeriods = nPeriods)
points(forecastObs[,,nPeriods],col="purple",pch=19,cex=0.5)
estForecastObs = sm.regression(forecastObs[,1,nPeriods],forecastObs[,2,nPeriods],display="none")
lines(estForecastObs$eval.points,estForecastObs$estimate,col="purple",lwd=3)

# forecast della est 
# forecastEst = forecastDiscrete(cbind(estFirst$eval.points,estFirst$estimate), est_field, speedFactor = speedFactor, nPeriods = nPeriods)
# lines(forecastEst[,1,nPeriods],forecastEst[,2,nPeriods],col="green",lwd=4)

legend("bottomright",legend=c("1960","2015","Forecast"),
       col=c("blue","red","purple"),lwd=c(3,3,3),bg="white")
dev.copy2pdf(file="Forecast19602015NonOverlapping.pdf",width=7,height=7,family = "mono")
par(op)

## plot forecast ----
lengthArrows = 1.0
latestYear = max(data$Year.end)
firstYear = min(data$Year.ini)
dataLatest = filter(data,Year.end == latestYear)
dataFirst = filter(data,Year.ini == firstYear)


dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 

# campo
signifVFest = significanceVF(est_field,X0,X1,0.05)
plot(est_field$x, type = "n", xlab = "GDP per capita (PPP in million 2017 USD, log scale)", ylab = "Life Expectancy", main = "",xaxt="n",xlim=c(6,12),ylim=c(50,90))
axis(1, at=6:12, labels=as.character(round(exp(6:12),-3)))
grid()

# arrows(est_field$x[,1], est_field$x[,2],
#        est_field$x[,1] + lengthArrows*est_field$estimator[,1], 
#        est_field$x[,2] + lengthArrows*est_field$estimator[,2],
#        length = 0.05, angle = 15, col = "black")
# punti
# points(dataFirst$GDP.t0,dataFirst$LE.t0,pch=19,col="blue",cex=0.5)
# text(dataFirst$GDP.t0,dataFirst$LE.t0,labels=dataFirst$countryCode,cex=1,pos=4,col="black")
# points(dataLatest$GDP.t0,dataLatest$LE.t0,pch=19,col="red",cex=0.5)
# text(dataLatest$GDP.t0,dataLatest$LE.t0,labels=dataLatest$countryCode,cex=0.5,pos=4,col="red")

# non parametriche
# estFirst = sm.regression(dataFirst$GDP.t0,dataFirst$LE.t0,display="none")
estLatest = sm.regression(dataLatest$GDP.t1,dataLatest$LE.t1,display="none")
# lines(estFirst$eval.points,estFirst$estimate,col="blue",lwd=3)
lines(estLatest$eval.points,estLatest$estimate,col="red",lwd=3)


# forecast delle obs 10
speedFactor = 0.1
nPeriods = 10/5 * 1/speedFactor
forecastObs = forecastDiscrete(cbind(dataLatest$GDP.t0,dataLatest$LE.t0), est_field, speedFactor = speedFactor, nPeriods = nPeriods)
# points(forecastObs[,,nPeriods],col="lightgreen",pch=19,cex=0.5)
estForecastObs = sm.regression(forecastObs[,1,nPeriods],forecastObs[,2,nPeriods],display="none")
lines(estForecastObs$eval.points,estForecastObs$estimate,col="lightgreen",lwd=3)

# forecast delle obs 20
speedFactor = 0.1
nPeriods = 20/5 * 1/speedFactor
forecastObs = forecastDiscrete(cbind(dataLatest$GDP.t0,dataLatest$LE.t0), est_field, speedFactor = speedFactor, nPeriods = nPeriods)
# points(forecastObs[,,nPeriods],col="forestgreen",pch=19,cex=0.5)
estForecastObs = sm.regression(forecastObs[,1,nPeriods],forecastObs[,2,nPeriods],display="none")
lines(estForecastObs$eval.points,estForecastObs$estimate,col="forestgreen",lwd=3)

# forecast delle obs 30
speedFactor = 0.1
nPeriods = 30/5 * 1/speedFactor
forecastObs = forecastDiscrete(cbind(dataLatest$GDP.t0,dataLatest$LE.t0), est_field, speedFactor = speedFactor, nPeriods = nPeriods)
points(forecastObs[,,nPeriods],col="darkgreen",pch=19,cex=0.5)
text(forecastObs[,,nPeriods],labels=dataLatest$countryCode,cex=0.75,pos=4,col="darkgreen")
estForecastObs = sm.regression(forecastObs[,1,nPeriods],forecastObs[,2,nPeriods],display="none")
lines(estForecastObs$eval.points,estForecastObs$estimate,col="darkgreen",lwd=3)

# forecast della est 
# forecastEst = forecastDiscrete(cbind(estFirst$eval.points,estFirst$estimate), est_field, speedFactor = speedFactor, nPeriods = nPeriods)
# lines(forecastEst[,1,nPeriods],forecastEst[,2,nPeriods],col="green",lwd=4)
legend("bottomright",legend=c("2015","2025 (Forecast)","2035 (Forecast)","2045 (Forecast)"),
       col=c("red","lightgreen","forestgreen","darkgreen"),lwd=c(3,3,3,3),bg="white")
dev.copy2pdf(file="Forecast19602015NonOverlapping102030.pdf",width=7,height=7,family = "mono")
par(op)


## plot direzioni verticali orizzontali campo vettoriale ----
library(fields)

dev.new()
signifVFest$signif[signifVFest$signif==0] = NA
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
Fx = signifVFest$signif*est_field$estimator[,1]/5
image.plot(x = unique(est_field$x[,1]), y = unique(est_field$x[,2]), z = matrix(Fx, nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="GDP per capita (PPP in million 2017 USD, log scale)",ylab="Life Expectancy",main="",xaxt="n")
axis(1, at=6:12, labels=as.character(round(exp(6:12),-3)))
points(data$GDP.t0,data$LE.t0,pch=19,col="black",cex=0.1)
dev.copy2pdf(file="Fx19602015NonOverlapping.pdf")

dev.new()
signifVFest$signif[signifVFest$signif==0] = NA
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
Fy = signifVFest$signif*est_field$estimator[,2]/5
image.plot(x = unique(est_field$x[,1]), y = unique(est_field$x[,2]), z = matrix(Fy, nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="GDP per capita (PPP in million 2017 USD, log scale)",ylab="Life Expectancy",main="",xaxt="n")
axis(1, at=6:12, labels=as.character(round(exp(6:12),-3)))
points(data$GDP.t0,data$LE.t0,pch=19,col="black",cex=0.1)
dev.copy2pdf(file="Fy19602015NonOverlapping.pdf")

dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
signifVFest$signif[signifVFest$signif==0] = NA
FyOverFx = signifVFest$signif*est_field$estimator[,2]/est_field$estimator[,1]
image.plot(x = unique(est_field$x[,1]), y = unique(est_field$x[,2]), z = matrix(log(FyOverFx), nrow=sqrt(nEval), ncol=sqrt(nEval)),xlab="GDP per capita (PPP in million 2017 USD, log scale)",ylab="Life Expectancy",main="",xaxt="n",nlevels=50)
axis(1, at=6:12, labels=as.character(round(exp(6:12),-3)))
points(data$GDP.t0,data$LE.t0,pch=19,col="black",cex=0.1)
dev.copy2pdf(file="FyoverFx19602015NonOverlapping.pdf")


# VAR a mano, ma forse giusto ----
# data = dplyr::filter(data, Year.ini >= 1980)
lm1 = lm(GDP.t1 ~ GDP.t0 + LE.t0, data = data)
lm2 = lm(LE.t1 ~ GDP.t0 + LE.t0, data = data)
A = matrix(c(coef(lm1)[2:3],coef(lm2)[2:3]),nrow=2,byrow=TRUE)
C = c(coef(lm1)[1],coef(lm2)[1])

# lm1 = lm(GDP.t1 ~ GDP.t0 + LE.t0 + 0, data = data)
# lm2 = lm(LE.t1 ~ GDP.t0 + LE.t0 + 0, data = data)
# A = matrix(c(coef(lm1)[1:2],coef(lm2)[1:2]),nrow=2,byrow=TRUE)
# C = c(0,0)

VFVAR = function(x){
  return((A-diag(2)) %*% x + C)
}
VFVARx = t(apply(x, 1, VFVAR))



## plot ----
latestYear = max(data$Year.end)
firstYear = min(data$Year.ini)
dataLatest = dplyr::filter(data, Year.end == latestYear)
dataFirst = dplyr::filter(data,Year.ini == firstYear)

dev.new()
# campo
lengthArrowsVAR = 1
plot(est_field$x, type = "n", xlab = "GDP per capita", ylab = "Life Expectancy", main = "Estimated VAR deg 1")
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
nPeriods = (2015-1960)/5 * 1/speedFactor
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

legend("bottomright",legend=c("1960","2015","Forecasted Observed","Forecasted Estimated"),
       col=c("black","red","purple","green"),pch=c(19,19,19,19),lwd=c(4,4,4,4),bg="white")
dev.copy2pdf(file="LogForecastVAR19602015NonOverlapping.pdf")

obs = cbind(dataLatest$GDP.t1,dataLatest$LE.t1)
diffVF = as.matrix(forecastObsVF-obs)
diffVAR = as.matrix(forecastObsVAR-obs)

errVF = apply(diffVF, MARGIN=1, FUN=function(x) sqrt(sum(x^2)))
errVAR = apply(diffVAR, MARGIN=1, FUN=function(x) sqrt(sum(x^2)))

summary(errVF)
summary(errVAR)

hist(errVF, 100)
hist(errVAR, 100, add=T, col="red")