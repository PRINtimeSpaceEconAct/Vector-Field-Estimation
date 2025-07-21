# Clear workspace and load dependencies
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
# setwd("/Users/davidefiaschi/Library/CloudStorage/OneDrive-UniversityofPisa/Cristiano\ Ricci\'s\ files\ -\ timeSpaceEvolutionEcAct/RVF/R\ code/Vector\ Field\ Estimation")

rm(list = ls())
DEBUG = TRUE
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
dev.copy2pdf(file="scatterArrowsPreston19601965.pdf")

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
#                                      sparse=FALSE, gc=TRUE)
est_field = NWfield(X0, X1, x=x, kernel.type="epa",h=1.0,
                                        chunk_size=1000,
                                        sparse=FALSE, gc=TRUE
                    #,method.h="sj"
                    )


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

# con i FE ----
timeInterval = 10 # 5 or 10 years
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
countriesDavide = read.csv(file = "datasets/countryCodeAnalysis.csv")

# Filter countries in Davide's file
data = data %>% filter(countryCode %in% countriesDavide[,2])

summary(log(data$GDP.t1/data$GDP.t0))


# count countries
length(unique(data$countryCode))

dataIni = data %>% select(c(GDP = GDP.t0,LE = LE.t0, Year = Year.ini,countryCode))
dataFin = data %>% select(c(GDP = GDP.t1,LE = LE.t1, Year = Year.end,countryCode)) %>% filter(Year == max(Year))



data = rbind(dataIni,dataFin) %>% arrange(countryCode)

data$GDP = log(data$GDP)

nObs = length(unique(data$countryCode))
nT = length(unique(data$Year))
nEval = 2500

# plot histograms of GDP growth


## qui tutto stima ----
# create array 
X = array(NA,dim=c(nObs, 2, nT))
for (year in unique(data$Year)) {
  X[,,(year - min(data$Year))/timeInterval+1] = cbind(data$GDP[data$Year == year], data$LE[data$Year == year])
}

# eval points
xGrid = seq(from=min(X[,1,]), to=max(X[,1,]), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X[,2,]), to=max(X[,2,]), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

Filtered = within_transform(X, FE = TRUE, TE = TRUE,
                            uniform_weights = TRUE, nEval_chunk = nEval,
                            x = x, kernel.type = "gauss",
                            method.h = "silverman", chunk_size = 512)


X0_raw = Filtered$X0_raw_unrolled
X0_star = Filtered$X0_star_unrolled
Y1_star = Filtered$Y1_star_unrolled
Y2_star = Filtered$Y2_star_unrolled
Y1 = Filtered$Y1_unrolled
Y2 = Filtered$Y2_unrolled


derivative_estimator_1 = compute_derivative_term(X0_raw, X0_star, x=x,
                                                 kernel.type="gauss", D=NULL, 
                                                 method.h="silverman", h=NULL, lambda=NULL, 
                                                 sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y1_star)

derivative_estimator_2 = compute_derivative_term(X0_raw, X0_star, x=x,
                                                 kernel.type="gauss", D=NULL,
                                                 method.h="silverman", h=NULL, lambda=NULL,
                                                 sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y2_star)

derivative_obs_1 = compute_derivative_term(X0_raw, X_star=X0_star[,,rep(1,nrow(X0_raw))], x=X0_raw,
                                           kernel.type="gauss", D=NULL, 
                                           method.h="silverman", h=NULL, lambda=NULL, 
                                           sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y1_star[,rep(1,nrow(Y1_star))])

derivative_obs_2 = compute_derivative_term(X0_raw, X_star=X0_star[,,rep(1,nrow(X0_raw))], x=X0_raw,
                                           kernel.type="gauss", D=NULL, 
                                           method.h="silverman", h=NULL, lambda=NULL, 
                                           sparse=FALSE, gc=FALSE, chunk_size=512, Y=Y2_star[,rep(1,nrow(Y2_star))])


meanPoint = apply(X0_raw, MARGIN = c(2), FUN = sum)/((nT-1)*nObs)
iBest = which.min(sqrt((x[,1]-meanPoint[1])^2 + (x[,2]-meanPoint[2])^2))


m10 = compute_m0(X_unrolled=X0_raw, Y_unrolled=Y1, beta=derivative_obs_1$estimator, x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
m20 = compute_m0(X_unrolled=X0_raw, Y_unrolled=Y2, beta=derivative_obs_2$estimator, x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])

VF_hat1 = compute_m(X0_raw, x, beta=derivative_estimator_1$estimator, m_0=m10, x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
VF_hat2 = compute_m(X0_raw, x, beta=derivative_estimator_2$estimator, m_0=m20, x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])

# Stitch together the two components of the vector field
VF_hat = cbind(VF_hat1, VF_hat2)

VF_hat_old <- VF_hat

# Call the new encapsulated function from panel.R
panel_vf_results <- estimate_panel_vf(X,
                                      x = x,
                                      nEval = nEval,
                                      FE = TRUE,
                                      TE = TRUE,
                                      uniform_weights = TRUE,
                                      kernel.type = "gauss",
                                      method.h = "silverman",
                                      chunk_size = 512,
                                      sparse = FALSE,
                                      gc = FALSE)

VF_hat_new <- panel_vf_results$VF_hat

# Compare the results
cat("\n\n--- Comparison between old and new method ---\n")
cat("Summary of differences (Old - New):\n")
print(summary(VF_hat_old - VF_hat_new))
cat("\nAre they numerically equal?\n")
print(all.equal(VF_hat_old, VF_hat_new))
cat("---------------------------------------------\n\n")


# rimuovi la media
# VF_hat = cbind(VF_hat1, VF_hat2) - c(mean(Y1),mean(Y2))

xGrid = seq(from=min(X[,1,]), to=max(X[,1,]), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X[,2,]), to=max(X[,2,]), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

# remove arrows farther than 0.1 from any point in X0_raw
for (i in 1:nEval){
  if (min(sqrt((x[i,1]-X0_raw[,1])^2 + (x[i,2]-X0_raw[,2])^2)) > 0.75) {
    VF_hat[i,] = c(0, 0)
  }
}

x[,1] = x[,1]/median(filter(data,Year==2010)$GDP)
x[,2] = x[,2]/median(filter(data,Year==2010)$LE)
VF_hat[,1] = VF_hat[,1]/median(filter(data,Year==2010)$GDP)
VF_hat[,2] = VF_hat[,2]/median(filter(data,Year==2010)$LE)

# lenVF = sqrt(VF_hat[,1]^2 + VF_hat[,2]^2)
# VF_hat[lenVF>1.5,] = c(0,0)
# ind.rm = c(which(x[,2] > x[,1]+0.125),which(x[,2] < x[,1]-0.25))
# VF_hat[ind.rm,] = c(0,0)

lengthArrows=(5/timeInterval)*1e-1*0.5
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
plot(x, type = "n", xlab = "Log GDP per capita (PPP in million 2017 USD, relative to 2010)", ylab = "Life Expectancy (relative to 2010)", main = "")
abline(h=1,lty=3)
abline(v=1,lty=3)
arrows(x[,1],x[,2],x[,1]+lengthArrows*VF_hat[,1],x[,2]+lengthArrows*VF_hat[,2],angle=15,col="black",length=0.05)
dataFirst = filter(data,Year == min(Year))
dataLatest = filter(data,Year == max(Year))
points((dataLatest$GDP)/median(filter(data,Year==2010)$GDP),dataLatest$LE/median(filter(data,Year==2010)$LE),pch=19,col="red",cex=0.5)
text((dataLatest$GDP)/median(filter(data,Year==2010)$GDP),dataLatest$LE/median(filter(data,Year==2010)$LE),labels=dataLatest$countryCode,cex=0.5,pos=4,col="red")
X0_raw_relative = X0_raw
X0_raw_relative[,1] = X0_raw_relative[,1]/median(filter(data,Year==2010)$GDP)
X0_raw_relative[,2] = X0_raw_relative[,2]/median(filter(data,Year==2010)$LE)
est.dens = sm.density(X0_raw_relative,display="none")
contour(est.dens$eval.points[,1], est.dens$eval.points[,2], est.dens$estimate,add=T,col="purple")

# grid()
dev.copy2pdf(file="outpics/VF_GDP_LE.pdf",width=7,height=7,family = "mono")

# reconstruct FE
VF_hat1 = compute_m(X0_raw, X0_raw, beta=derivative_obs_1$estimator, m_0=m10, x0=x[iBest,], beta_0=derivative_estimator_1$estimator[iBest,])
VF_hat2 = compute_m(X0_raw, X0_raw, beta=derivative_obs_2$estimator, m_0=m20, x0=x[iBest,], beta_0=derivative_estimator_2$estimator[iBest,])

YObs = aperm(array(cbind(Y1,Y2),dim = c(nT-1,nObs,2)),c(2,3,1))
VFObs = aperm(array(cbind(VF_hat1,VF_hat2),dim = c(nT-1,nObs,2)),c(2,3,1))

alpha_i_hat = apply(YObs - VFObs, MARGIN = c(1, 2), FUN = sum)/(nT-1)
gamma_t_hat = t(apply(YObs - VFObs, MARGIN = c(2, 3), FUN = sum)/nObs)


# filtro FE e TE e ristimo
Y1_rolled = array(Y1,dim=c(nT-1,nObs))
Y2_rolled = array(Y2,dim=c(nT-1,nObs))
FE1Rep = t(array(rep(alpha_i_hat[,1],nT-1),dim=c(nObs,nT-1)))
FE2Rep = t(array(rep(alpha_i_hat[,2],nT-1),dim=c(nObs,nT-1)))
TE1Rep = array(rep(gamma_t_hat[,1],nObs),dim=c(nT-1,nObs))
TE2Rep = array(rep(gamma_t_hat[,2],nObs),dim=c(nT-1,nObs))

# FE1Rep = (array(rep(alpha_i_hat[,1],nT-1),dim=c(nT-1,nObs)))
# FE2Rep = (array(rep(alpha_i_hat[,2],nT-1),dim=c(nT-1,nObs)))
# TE1Rep = t(array(rep(gamma_t_hat[,1],nObs),dim=c(nObs,nT-1)))
# TE2Rep = t(array(rep(gamma_t_hat[,2],nObs),dim=c(nObs,nT-1)))


Y1_rolled_filtered = Y1_rolled - FE1Rep - TE1Rep
Y2_rolled_filtered = Y2_rolled - FE2Rep - TE2Rep

Y1_filtered = t(Y1_rolled_filtered)
dim(Y1_filtered) = c((nT-1) * nObs)
Y2_filtered = t(Y2_rolled_filtered)
dim(Y2_filtered) = c((nT-1) * nObs)

Y_filtered = cbind(Y1_filtered,Y2_filtered) 
X0 = X0_raw
X1 = Y_filtered + X0

DEBUG = FALSE
est_field_adaptive = NWfieldAdaptive(X0, X1, x=x, kernel.type="gauss",
                                     chunk_size=1000,
                                     sparse=FALSE, gc=TRUE, 
                                     hOpt = TRUE, alphaOpt = TRUE)
est_field_adaptive_bootstrap = bootstrapKernelFieldErrors(est_field_adaptive, B = 100)
signifBoot = significanceBootstrap(est_field_adaptive,est_field_adaptive_bootstrap)

# signifInd = which(signifBoot)
signifInd = 1:length(signifBoot)
# dev.new()
lengthArrows=0.1
plot(x, type = "n", xlab = "GDP (log)", ylab="LE", main = " ")
arrows(est_field_adaptive$x[signifInd,1], est_field_adaptive$x[signifInd,2],
       est_field_adaptive$x[signifInd,1] + lengthArrows*est_field_adaptive$estimator[signifInd,1], 
       est_field_adaptive$x[signifInd,2] + lengthArrows*est_field_adaptive$estimator[signifInd,2],
       length = 0.05, angle = 15, col = "black")
dataFirst = filter(data,Year == min(Year))
dataLatest = filter(data,Year == max(Year))
points((dataFirst$GDP),dataFirst$LE,pch=19,col="blue",cex=0.5)
points((dataLatest$GDP),dataLatest$LE,pch=19,col="red",cex=0.5)
text((dataLatest$GDP),dataLatest$LE,labels=dataLatest$countryCode,cex=0.5,pos=4,col="red")


## plot FE ----
GDP0 = arrange(filter(data,Year == min(Year)),countryCode)$GDP
LE0 = arrange(filter(data,Year == min(Year)),countryCode)$LE
names0 = arrange(filter(data,Year == min(Year)),countryCode)$countryCode

dev.new()
plot(GDP0,alpha_i_hat[,1], main = "", xlab = "GDP (log)", ylab = "FE", pch=19, cex = 0.5)
abline(h=0)
grid()
text(GDP0,alpha_i_hat[,1],labels=names0,cex=0.5,pos=4,col="red")
dev.copy2pdf(file="outpics/FE_GDP0.pdf",width=7,height=7,family = "mono")

dev.new()
plot(LE0,alpha_i_hat[,2], main = "", xlab = "LE", ylab = "FE", pch=19, cex = 0.5)
abline(h=0)
grid()
text(LE0,alpha_i_hat[,2],labels=names0,cex=0.5,pos=4,col="red")
dev.copy2pdf(file="outpics/FE_LE0.pdf",width=7,height=7,family = "mono")

## plot TE ----
times = sort(unique(data$Year))

dev.new()
plot(times[2:nT],gamma_t_hat[,1], main = "", xlab = "Years", ylab = "TE GDP (log)", pch=19, cex = 0.5,type="b")
abline(h=0)
grid()
dev.copy2pdf(file="outpics/TE_GDP0.pdf",width=7,height=7,family = "mono")

dev.new()
plot(times[2:nT],gamma_t_hat[,2], main = "", xlab = "Years", ylab = "TE LE", pch=19, cex = 0.5,type="b")
abline(h=0)
grid()
dev.copy2pdf(file="outpics/TE_LE0.pdf",width=7,height=7,family = "mono")

## plot VF_obs ----
dev.new()
lengthArrowsVF_obs = 0.1

# Add fixed effects to the vector field observations
VFObs_with_FE = sweep(VFObs, c(1, 2), alpha_i_hat, "+")
# Add time effects
VFObs_with_FE_TE = sweep(VFObs_with_FE, c(3, 2), gamma_t_hat, "+")
# Reshape to match X0_raw for plotting
VF_FE_TE_obs_unrolled = matrix(aperm(VFObs_with_FE_TE, c(3,1,2)), ncol=2)

plot(X0_raw, type = "n", xlab = "GDP (log)", ylab="LE", main = "Reconstructed VF at observed points")
arrows(X0_raw[,1], X0_raw[,2], X0_raw[,1] + lengthArrowsVF_obs*VF_FE_TE_obs_unrolled[,1], X0_raw[,2] + lengthArrowsVF_obs*VF_FE_TE_obs_unrolled[,2], angle=15, col="black", length=0.05)
dataFirst = filter(data,Year == min(Year))
dataLatest = filter(data,Year == max(Year))
points(dataFirst$GDP, dataFirst$LE, pch=19, col="blue", cex=0.5)
points(dataLatest$GDP, dataLatest$LE, pch=19, col="red", cex=0.5)
text(dataLatest$GDP, dataLatest$LE, labels=dataLatest$countryCode, cex=0.5, pos=4, col="red")
grid()
dev.copy2pdf(file="outpics/VF_reconstructed_obs_points.pdf", width=7, height=7, family="mono")

# Error plot ----
dev.new()
lengthArrowsError = 0.3
Error =YObs #- VFObs_with_FE_TE
Error_unrolled = matrix(aperm(Error, c(3,1,2)), ncol=2)
plot(X0_raw, type = "n", xlab = "GDP (log)", ylab="LE", main = "Error")
arrows(X0_raw[,1], X0_raw[,2], X0_raw[,1] + lengthArrowsError*Error_unrolled[,1], X0_raw[,2] + lengthArrowsError*Error_unrolled[,2], angle=15, col="black", length=0.05)
dataFirst = filter(data,Year == min(Year))
dataLatest = filter(data,Year == max(Year))
points(dataFirst$GDP, dataFirst$LE, pch=19, col="blue", cex=0.5)
points(dataLatest$GDP, dataLatest$LE, pch=19, col="red", cex=0.5)
text(dataLatest$GDP, dataLatest$LE, labels=dataLatest$countryCode, cex=0.5, pos=4, col="red")
grid()

# analisi dati ----
data = dataset.GDP.LE
data = data %>% mutate(LogGDP.t0 = log(GDP.t0), LogGDP.t1 = log(GDP.t1))
data = data %>% mutate(GDPGR = LogGDP.t1-LogGDP.t0)
dataOutlier = filter(data, GDPGR > 0.4 | GDPGR < -0.4)
plot(dataOutlier$Year.ini, dataOutlier$GDPGR,pch=19,cex=0.5,col="black",xlab="Year",ylab="GDP GR (E` il cumulato a 5 anni)",main="")
text(dataOutlier$Year.ini, dataOutlier$GDPGR,labels = dataOutlier$countryCode,cex=0.5,pos=4,col="red")


