rm(list=ls())

setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
# setwd("/Users/davidefiaschi/Library/CloudStorage/OneDrive-UniversityofPisa/Cristiano\ Ricci\'s\ files\ -\ timeSpaceEvolutionEcAct/RVF/R\ code/Vector\ Field\ Estimation")

rm(list = ls())
DEBUG = FALSE
source("src/libs/loadLib.R")
library(dplyr)
library(sm)
library(xtable)
library(latex2exp)

load("datasets/datasetGDP.LE.NonOverlapping.RData")
data = dataset.GDP.LE
countriesDavide = read.csv(file = "../Vector Field Estimation/datasets/countryCodeAnalysis.csv")

# Filter countries in Davide's file
data = data %>% filter(countryCode %in% countriesDavide[,2])

# timeInterval 10 
timeInterval = 5
# df1 <- data
# df2 <- data

# data_10y <- inner_join(
#     df1, 
#     df2, 
#     by = "countryCode",
#     suffix = c(".x", ".y"),
#     relationship = "many-to-many"
# ) %>%
#     filter(Year.end.x == Year.ini.y) %>%
#     transmute(
#         countryCode = countryCode,
#         Year.ini = Year.ini.x,
#         Year.end = Year.end.y,
#         GDP.t0 = GDP.t0.x,
#         LE.t0 = LE.t0.x,
#         GDP.t1 = GDP.t1.y,
#         LE.t1 = LE.t1.y
#     )
# 
# # Filter to get non-overlapping 10-year intervals (e.g., 1960-1970, 1970-1980)
# startYearOfDecade = min(data_10y$Year.ini) %% 10
# data = data_10y %>% filter(Year.ini %% 10 == startYearOfDecade)

data = data %>% mutate(GDP = log(GDP.t0), Delta = (log(GDP.t1)-log(GDP.t0))/5) %>% 
    select(countryCode, Year = Year.ini, GDP, Delta)


data0 = filter(data, Year < max(data$Year)) %>% arrange(Year,countryCode)
data1 = filter(data, Year > min(data$Year)) %>% arrange(Year,countryCode)
dataYearIni = filter(data, Year == min(data$Year))
dataYearEnd = filter(data, Year == max(data$Year))

# regressione 
plot(data$GDP,data$Delta,xlab="logGDP",ylab="Delta logGDP 10Anni")
sm.regression(data$GDP,data$Delta)

X0 = cbind(data0$GDP,data0$Delta)
X1 = cbind(data1$GDP,data1$Delta)

# parameters ----
nEval = 2500

xGrid = seq(from=0.9*min(c(X0[,1],X1[,1])), to=1.1*max(c(X0[,1],X1[,1])), length.out=round(sqrt(nEval)))
yGrid = seq(from=0.9*min(c(X0[,2],X1[,2])), to=1.1*max(c(X0[,2],X1[,2])), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))


est_field_adaptive = NWfieldAdaptive(X0, X1, x=x, kernel.type="epa",method.h = "silverman",
                                     chunk_size=1000,
                                     gc=TRUE,
                                     hOpt = TRUE, alphaOpt = TRUE)
# est_field_adaptive = NWfieldAdaptive(X0, X1, x=x, kernel.type="epa",
#                                      chunk_size=1000,h=0.1,
#                                      gc=TRUE, 
#                                      hOpt = FALSE, alphaOpt = FALSE)

# bootstrap
est_field_adaptive_bootstrap = bootstrapKernelFieldErrors(est_field_adaptive, B = 100)
signifBoot = significanceBootstrap(est_field_adaptive,est_field_adaptive_bootstrap)


## plot campo stimato ----
dev.new()
op <- par(family = "mono") 
signifVFest = significanceVF(est_field_adaptive)
lengthArrows=0.1
plot(x, type = "n", xlab="logGDP",ylab="Delta logGDP 10 Anni", main = "Tutte frecce",ylim=c(-0.1,0.1))
points(dataYearIni$GDP,dataYearIni$Delta,col="blue")
points(dataYearEnd$GDP,dataYearEnd$Delta,col="red")
sm.reg = sm.regression(data$GDP,data$Delta,display="none")
lines(sm.reg$eval.points,sm.reg$estimate,col="purple")
arrows(est_field_adaptive$x[,1], est_field_adaptive$x[,2],
       est_field_adaptive$x[,1] + lengthArrows*est_field_adaptive$estimator[,1], 
       est_field_adaptive$x[,2] + lengthArrows*est_field_adaptive$estimator[,2],
       length = 0.05, angle = 15, col = "black")
abline(h=0)
abline(v=0)
par(op)
dev.copy2pdf(file="fiaschilavezzi.pdf")

## plot campo stimato significativo ----
signifInd = which(signifBoot)
dev.new()
op <- par(family = "mono")
lengthArrows=0.2
plot(x, type = "n", xlab="logGDP",ylab="Delta logGDP 10 Anni", main = "Solo significative",ylim=c(-0.05,0.075))
points(dataYearIni$GDP,dataYearIni$Delta,col="blue",pch=19,cex=.45)
text(dataYearIni$GDP,dataYearIni$Delta,dataYearIni$countryCode,col="blue",pos=4,cex=.45)
sm.reg = sm.regression(data$GDP,data$Delta,display="none")
lines(sm.reg$eval.points,sm.reg$estimate,col="purple")
arrows(est_field_adaptive$x[signifInd,1], est_field_adaptive$x[signifInd,2],
       est_field_adaptive$x[signifInd,1] + lengthArrows*est_field_adaptive$estimator[signifInd,1],
       est_field_adaptive$x[signifInd,2] + lengthArrows*est_field_adaptive$estimator[signifInd,2],
       length = 0.05, angle = 15, col = "black")
abline(h=mean(data$Delta))
abline(h=0)
abline(v=0)
par(op)
dev.copy2pdf(file="fiaschilavezziSignif.pdf")
# 
# 
