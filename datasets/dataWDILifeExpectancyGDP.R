#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Dataset of GDP per capita versus life expectancy at country level
#
# Source: Wordl Bank
#
# Update: February 13, 2025
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rm(list = ls())

#Libraries ####
library(WDI)
library(pwt10)
library(sm)
#library(KernSmooth)
#library(locfit)
library(locpol)

data("pwt10.01")

# Download life expectancy at birth from World Bank dataset ####
# (https://data.worldbank.org/indicator/)

#List of variables available at the World Bank dataset
#WDIsearch()

## Life expectancy at birth ####
life.exp.at.birth <- WDI(country = "all",indicator="SP.DYN.LE00.IN",start=1960,end = NULL,extra=TRUE)

#Eliminate the aggregates among countries
iii <- !is.na(life.exp.at.birth$capital) & !(life.exp.at.birth$capital=="")
life.exp.at.birth <- life.exp.at.birth[iii,]

## GDP per capita, PPP (constant 2021 international $) (https://data.worldbank.org/indicator/NY.GDP.PCAP.PP.KD) ####
#For Penn World Table we use: "Expenditure-side real GDP at chained PPPs (in million 2017 USD)."
GDP.per.capita.PPP <- WDI(country = "all",indicator="NY.GDP.PCAP.PP.KD",start=1960,end = NULL,extra=TRUE)

iii <- !is.na(GDP.per.capita.PPP$capital) & !(GDP.per.capita.PPP$capital=="")
GDP.per.capita.PPP <- GDP.per.capita.PPP[iii,]

#List of countries ####
#It the intersection between countries in World Bank and Penn World Table (the last one is the longest dataset on GDP per capita)
list.countries <- sort(intersect(unique(life.exp.at.birth$iso3c),unique(pwt10.01$isocode)))

#List of years ####
list.years <- sort(unique(life.exp.at.birth$year))
#Pick up 2023 because no data for life expectancy
list.years <- list.years[-length(list.years)]

#Bulding the dataset as a balanced panel of countries ####
#Initializing the matrix for life expecntancy and GDP per capita
life.exp.matrix <- matrix(NA,nrow=length(list.countries),ncol=length(list.years))
GDP.per.capita.matrix <- matrix(NA,nrow=length(list.countries),ncol=length(list.years))
GDP.per.capita.PWT.matrix <- matrix(NA,nrow=length(list.countries),ncol=length(list.years))

colnames(life.exp.matrix) <- list.years
colnames(GDP.per.capita.matrix) <- list.years
colnames(GDP.per.capita.PWT.matrix) <- list.years
rownames(life.exp.matrix) <- list.countries
rownames(GDP.per.capita.matrix) <- list.countries
rownames(GDP.per.capita.PWT.matrix) <- list.countries

for (i in 1:length(list.countries)){
  
  iii.LE <- life.exp.at.birth$iso3c==list.countries[i]
  iii.GDP <- GDP.per.capita.PPP$iso3c==list.countries[i]
  iii.GDP.PWT <- pwt10.01$isocode==list.countries[i]
  
  for (j in 1:length(list.years)){
    
    jjj.LE <- life.exp.at.birth$year==list.years[j]
    jjj.GDP <- GDP.per.capita.PPP$year==list.years[j]
    jjj.GDP.PWT <- pwt10.01$year==list.years[j]
    
    life.exp.matrix[i,j] <- life.exp.at.birth$SP.DYN.LE00.IN[iii.LE & jjj.LE]
    GDP.per.capita.matrix[i,j] <- GDP.per.capita.PPP$NY.GDP.PCAP.PP.KD[iii.GDP & jjj.GDP]
    if (sum(iii.GDP.PWT & jjj.GDP.PWT)==0){}else{
      GDP.per.capita.PWT.matrix[i,j] <- pwt10.01$rgdpe[iii.GDP.PWT & jjj.GDP.PWT]/pwt10.01$pop[iii.GDP.PWT & jjj.GDP.PWT] 
    }
  }
}

## Merging World Bank and PWT data (before 2019 PWT, after WB. NB: the base years is different, so they are not directly usable) ####

#summary(GDP.per.capita.PWT.matrix[,which(list.years==1990)]/GDP.per.capita.matrix[,which(list.years==1990)])
GDP.per.capita.merge <- cbind(GDP.per.capita.PWT.matrix[,1:which(list.years==2019)],GDP.per.capita.matrix[,(which(list.years==2019)+1):ncol(GDP.per.capita.matrix)])

#Selection of the first year of the balenced dataset (there exists a strong trade-off before and after 1970; 2019 is the last year in Penn World Table)
starting.year <- 1960
final.year <- 2019

iii.sel.GDP <- !is.na(rowSums(GDP.per.capita.merge[,(which(list.years==starting.year)):(which(list.years==final.year))]))

iii.sel.LE <- !is.na(rowSums(life.exp.matrix[,(which(list.years==starting.year)):(which(list.years==final.year))]))

iii.sel <- iii.sel.GDP & iii.sel.LE

#Pick up Ireland and Luxembourg as tax heaven
iii.sel[which(list.countries=="IRL")] <- FALSE
iii.sel[which(list.countries=="LUX")] <- FALSE

GDP.per.capita.sel <- log(GDP.per.capita.merge[iii.sel,(which(list.years==starting.year)):(which(list.years==final.year))])
life.exp.sel <- life.exp.matrix[iii.sel,(which(list.years==starting.year)):(which(list.years==final.year))]

est.1960 <- sm.regression(GDP.per.capita.sel[,1],life.exp.sel[,1],display="none")
est.1970 <- sm.regression(GDP.per.capita.sel[,11],life.exp.sel[,11],display="none")
est.1980 <- sm.regression(GDP.per.capita.sel[,21],life.exp.sel[,21],display="none")
est.1990 <- sm.regression(GDP.per.capita.sel[,31],life.exp.sel[,31],display="none")
est.2000 <- sm.regression(GDP.per.capita.sel[,41],life.exp.sel[,41],display="none")
est.2010 <- sm.regression(GDP.per.capita.sel[,51],life.exp.sel[,51],display="none")
est.2015 <- sm.regression(GDP.per.capita.sel[,(ncol(GDP.per.capita.sel)-4)],life.exp.sel[,ncol(life.exp.sel)],display="none")

#h <- dpill(GDP.per.capita.sel[,1], life.exp.sel[,1])
#fitted.1960 <- locpoly(x= GDP.per.capita.sel[,1], y = life.exp.sel[,1],degree=1,bandwidth = h)
#lines(fitted.1960$x,fitted.1960$y,col="red",lwd=2)
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 

plot(GDP.per.capita.sel[,1],life.exp.sel[,1],pch=19,cex=0.25,col="blue",ylim=c(30,85),xlim=c(6.,11.5),ylab="Life expectancy at birth",xlab="GDP per capita (PPP in million 2017 USD, log scale)")
points(GDP.per.capita.sel[,(ncol(GDP.per.capita.sel)-4)],life.exp.sel[,ncol(life.exp.sel)],pch=19,cex=0.25,col="red")

lines(est.2015$eval.points,est.2015$estimate,col="red",lwd=2)
lines(est.1960$eval.points,est.1960$estimate,col="blue",lwd=2)
lines(est.1970$eval.points,est.1970$estimate,col="gray40",lwd=1.5)
lines(est.1980$eval.points,est.1980$estimate,col="gray30",lwd=1.5)
lines(est.1990$eval.points,est.1990$estimate,col="gray20",lwd=1.5)
lines(est.2000$eval.points,est.2000$estimate,col="gray10",lwd=1.5)
lines(est.2010$eval.points,est.2010$estimate,col="gray0",lwd=1.5)
grid()
legend("bottomright",c("1960","1970","1980","1990","2000","2010","2015"),col=c("blue","gray40","gray30","gray20","gray10","gray0","red"),lty=1,lwd=1.5)

dev.copy2eps(file="prestonCurveDifferentYears.eps",width=7,height=7,family = "mono")
par(op)


#Relative values ####
GDP.per.capita.sel.REL <- GDP.per.capita.sel/matrix(colMeans(GDP.per.capita.sel),ncol=ncol(GDP.per.capita.sel),nrow=nrow(GDP.per.capita.sel),byrow=T)
life.exp.sel.REL <- life.exp.sel/matrix(colMeans(life.exp.sel),ncol=ncol(life.exp.sel),nrow=nrow(life.exp.sel),byrow=T)


est.1960 <- sm.regression(GDP.per.capita.sel.REL[,1],life.exp.sel.REL[,1],display="none")
est.1970 <- sm.regression(GDP.per.capita.sel.REL[,11],life.exp.sel.REL[,11],display="none")
est.1980 <- sm.regression(GDP.per.capita.sel.REL[,21],life.exp.sel.REL[,21],display="none")
est.1990 <- sm.regression(GDP.per.capita.sel.REL[,31],life.exp.sel.REL[,31],display="none")
est.2000 <- sm.regression(GDP.per.capita.sel.REL[,41],life.exp.sel.REL[,41],display="none")
est.2010 <- sm.regression(GDP.per.capita.sel.REL[,51],life.exp.sel.REL[,51],display="none")
est.2019 <- sm.regression(GDP.per.capita.sel.REL[,ncol(GDP.per.capita.sel.REL)],life.exp.sel.REL[,ncol(life.exp.sel.REL)],display="none")

#h <- dpill(GDP.per.capita.sel.REL[,1], life.exp.sel.REL[,1])
#fitted.1960 <- locpoly(x= GDP.per.capita.sel.REL[,1], y = life.exp.sel.REL[,1],degree=1,bandwidth = h)
#lines(fitted.1960$x,fitted.1960$y,col="red",lwd=2)

plot(GDP.per.capita.sel.REL[,1],life.exp.sel.REL[,1],pch=19,cex=0.25,col="red",ylim=c(0.2,1.4),xlim=c(0,4),ylab="Life expectancy",xlab="GDP per capita PPP (international $)")
points(GDP.per.capita.sel.REL[,ncol(GDP.per.capita.sel.REL)],life.exp.sel.REL[,ncol(life.exp.sel.REL)],pch=19,cex=0.25)

lines(est.2019$eval.points,est.2019$estimate,col="black",lwd=2)
lines(est.1960$eval.points,est.1960$estimate,col="red",lwd=2)
lines(est.1970$eval.points,est.1970$estimate,col="lightgray",lwd=1)
lines(est.1980$eval.points,est.1980$estimate,col="gray",lwd=1.5)
lines(est.1990$eval.points,est.1990$estimate,col="gray",lwd=2)
lines(est.2000$eval.points,est.2000$estimate,col="darkgray",lwd=2.5)
#lines(est.2010$eval.points,est.2010$estimate,col="gray",lwd=3)
grid()

save(GDP.per.capita.sel.REL,life.exp.sel.REL,life.exp.sel,GDP.per.capita.sel,file="lifeExpGDPPerCapita19602019.RData")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Building a dataset with 5 year lag  ----
library(dplyr)
load("lifeExpGDPPerCapita19602019.RData")

# to remove years after 2015 in order to avoid a transition of 4 years instead of 5.
GDP.per.capita.sel = GDP.per.capita.sel[,colnames(GDP.per.capita.sel) <= 2015]
GDP.per.capita.sel.REL = GDP.per.capita.sel.REL[,colnames(GDP.per.capita.sel.REL) <= 2015]
life.exp.sel = life.exp.sel[,colnames(life.exp.sel) <= 2015]
life.exp.sel.REL = life.exp.sel.REL[,colnames(life.exp.sel.REL) <= 2015]
    
starting.year = as.numeric(colnames(GDP.per.capita.sel.REL)[1])
final.year = as.numeric(colnames(GDP.per.capita.sel.REL)[ncol(GDP.per.capita.sel.REL)])

#The length of lag
lag.length <- 5

## The dataset with not overlapping observations ----
yearsIntervals <- seq(from=1,to=ncol(GDP.per.capita.sel.REL),by=lag.length)
## The last year is imposed in the list of years, in case is not mutiple of lag.length
# yearsIntervals[length(yearsIntervals)+1] <- ncol(GDP.per.capita.sel.REL)

#GDP per capita
GDP.per.capita.t0 <- GDP.per.capita.sel[,yearsIntervals[-length(yearsIntervals)]]
GDP.per.capita.t1 <- GDP.per.capita.sel[,yearsIntervals[-1]]
#Life expectancy at birth
life.exp.t0 <- life.exp.sel[,yearsIntervals[-length(yearsIntervals)]]
life.exp.t1 <- life.exp.sel[,yearsIntervals[-1]]

#Relative GDP per capita
GDP.per.capita.REL.t0 <- GDP.per.capita.sel.REL[,yearsIntervals[-length(yearsIntervals)]]
GDP.per.capita.REL.t1 <- GDP.per.capita.sel.REL[,yearsIntervals[-1]]
#Relative life expectancy at birth
life.exp.REL.t0 <- life.exp.sel.REL[,yearsIntervals[-length(yearsIntervals)]]
life.exp.REL.t1 <- life.exp.sel.REL[,yearsIntervals[-1]]
year.ini <- (starting.year:final.year)[yearsIntervals[-length(yearsIntervals)]]
year.end <- year.ini + lag.length

# dataset in formmat RVF friendly 
dataset.GDP.LE <- cbind(c(GDP.per.capita.t0),c(life.exp.t0),c(GDP.per.capita.t1),c(life.exp.t1),c(matrix(year.ini,nrow = nrow(GDP.per.capita.t0),ncol=ncol(GDP.per.capita.t0),byrow=T)),c(matrix(year.end,nrow = nrow(GDP.per.capita.t0),ncol=ncol(GDP.per.capita.t0),byrow=T)))
colnames(dataset.GDP.LE) <- c("GDP.t0","LE.t0","GDP.t1","LE.t1","Year.ini","Year.end")
rownames(dataset.GDP.LE) <- rep(rownames(GDP.per.capita.sel),times=(length(yearsIntervals)-1))

dataset.GDP.LE = data.frame(dataset.GDP.LE)
dataset.GDP.LE = dataset.GDP.LE %>% mutate(countryCode = rownames(dataset.GDP.LE)) %>% mutate(countryCode = substr(countryCode, 1, 3))
rownames(dataset.GDP.LE) = NULL


# dataset in formmat RVF friendly, REL
dataset.GDP.LE.REL <- cbind(c(GDP.per.capita.REL.t0),c(life.exp.REL.t0),c(GDP.per.capita.REL.t1),c(life.exp.REL.t1),c(matrix(year.ini,nrow = nrow(GDP.per.capita.REL.t0),ncol=ncol(GDP.per.capita.REL.t0),byrow=T)),c(matrix(year.end,nrow = nrow(GDP.per.capita.REL.t0),ncol=ncol(GDP.per.capita.REL.t0),byrow=T)))
colnames(dataset.GDP.LE.REL) <- c("GDP.REL.t0","LE.REL.t0","GDP.REL.t1","LE.REL.t1","Year.ini","Year.end")
rownames(dataset.GDP.LE.REL) <- rep(rownames(GDP.per.capita.sel.REL),times=(length(yearsIntervals)-1))

dataset.GDP.LE.REL = data.frame(dataset.GDP.LE.REL)
dataset.GDP.LE.REL = dataset.GDP.LE.REL %>% mutate(countryCode = rownames(dataset.GDP.LE.REL)) %>% mutate(countryCode = substr(countryCode, 1, 3))
rownames(dataset.GDP.LE.REL) = NULL

save(dataset.GDP.LE,dataset.GDP.LE.REL,file="datasetGDP.LE.NonOverlapping.RData")

## The dataset with overlapping observations ----
yearsIndexInitial <- seq(from=1,to=ncol(GDP.per.capita.sel.REL)-lag.length,by=1)
yearsIndexEnd = yearsIndexInitial + lag.length

GDP.per.capita.t0 <- GDP.per.capita.sel[,yearsIndexInitial]
GDP.per.capita.t1 <- GDP.per.capita.sel[,yearsIndexEnd]
#Life expectancy at birth
life.exp.t0 <- life.exp.sel[,yearsIndexInitial]
life.exp.t1 <- life.exp.sel[,yearsIndexEnd]

#Relative GDP per capita
GDP.per.capita.REL.t0 <- GDP.per.capita.sel.REL[,yearsIndexInitial]
GDP.per.capita.REL.t1 <- GDP.per.capita.sel.REL[,yearsIndexEnd]
#Relative life expectancy at birth
life.exp.REL.t0 <- life.exp.sel.REL[,yearsIndexInitial]
life.exp.REL.t1 <- life.exp.sel.REL[,yearsIndexEnd]
year.ini <- colnames(GDP.per.capita.sel)[yearsIndexInitial]
year.end <- colnames(GDP.per.capita.sel)[yearsIndexEnd]

# dataset in formmat RVF friendly 
dataset.GDP.LE <- cbind(c(GDP.per.capita.t0),c(life.exp.t0),c(GDP.per.capita.t1),c(life.exp.t1),c(matrix(year.ini,nrow = nrow(GDP.per.capita.t0),ncol=ncol(GDP.per.capita.t0),byrow=T)),c(matrix(year.end,nrow = nrow(GDP.per.capita.t0),ncol=ncol(GDP.per.capita.t0),byrow=T)))
colnames(dataset.GDP.LE) <- c("GDP.t0","LE.t0","GDP.t1","LE.t1","Year.ini","Year.end")
rownames(dataset.GDP.LE) <- rep(rownames(GDP.per.capita.sel),times=(length(yearsIndexInitial)))

dataset.GDP.LE = data.frame(dataset.GDP.LE)
dataset.GDP.LE = dataset.GDP.LE %>% mutate(countryCode = rownames(dataset.GDP.LE)) %>% mutate(countryCode = substr(countryCode, 1, 3))
rownames(dataset.GDP.LE) = NULL

# dataset in formmat RVF friendly, REL
dataset.GDP.LE.REL <- cbind(c(GDP.per.capita.REL.t0),c(life.exp.REL.t0),c(GDP.per.capita.REL.t1),c(life.exp.REL.t1),c(matrix(year.ini,nrow = nrow(GDP.per.capita.REL.t0),ncol=ncol(GDP.per.capita.REL.t0),byrow=T)),c(matrix(year.end,nrow = nrow(GDP.per.capita.t0),ncol=ncol(GDP.per.capita.t0),byrow=T)))
colnames(dataset.GDP.LE.REL) <- c("GDP.REL.t0","LE.REL.t0","GDP.REL.t1","LE.REL.t1","Year.ini","Year.end")
rownames(dataset.GDP.LE.REL) <- rep(rownames(GDP.per.capita.sel.REL),times=(length(yearsIndexInitial)))

dataset.GDP.LE.REL = data.frame(dataset.GDP.LE.REL)
dataset.GDP.LE.REL = dataset.GDP.LE.REL %>% mutate(countryCode = rownames(dataset.GDP.LE.REL)) %>% mutate(countryCode = substr(countryCode, 1, 3))
rownames(dataset.GDP.LE.REL) = NULL

save(dataset.GDP.LE,dataset.GDP.LE.REL,file="datasetGDP.LE.Overlapping.RData")

