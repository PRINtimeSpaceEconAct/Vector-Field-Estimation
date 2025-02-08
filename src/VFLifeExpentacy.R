# Clear workspace and load dependencies
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
rm(list = ls())
DEBUG = TRUE
source("src/libs/loadLib.R")

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
xGrid = seq(from=min(c(X0[,1],X1[,1])), to=max(c(X0[,1],X1[,1])), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(c(X0[,2],X1[,2])), to=max(c(X0[,2],X1[,2])), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))



# est_field = LLfieldAdaptive(X0, X1, x=x, kernel.type="epa",method.h = "sj",
#                                      chunk_size=1000,
#                                      sparse=FALSE, gc=TRUE, alpha=0.5)

est_field = LLfield(X0, X1, x=x, kernel.type="epa",method.h = "sj",
                                     chunk_size=1000,
                                     sparse=FALSE, gc=TRUE)


# est_fieldTREND = NWfield(X0, X1, x=x, kernel.type="gauss",h = 1000000,
#                     chunk_size=1000,
#                     sparse=FALSE, gc=TRUE)
# 
# est_field$estimator = est_field$estimator - est_fieldTREND$estimator

## plot campo stimato ----
lengthArrows = 1.0
latestYear = max(data$Year.end)
dataLatest = filter(data, Year.end == latestYear)

# dev.new()
plot(est_field$x, type = "n", xlab = "GDP per capita", ylab = "Life Expectancy", main = "Estimated Vector Field")
arrows(est_field$x[,1], est_field$x[,2],
       est_field$x[,1] + lengthArrows*est_field$estimator[,1], 
       est_field$x[,2] + lengthArrows*est_field$estimator[,2],
       length = 0.05, angle = 15, col = "blue")
points(dataLatest$GDP.t1,dataLatest$LE.t1)
text(dataLatest$GDP.t1,dataLatest$LE.t1,labels=dataLatest$countryCode, col="red",cex=0.5,pos=4)




