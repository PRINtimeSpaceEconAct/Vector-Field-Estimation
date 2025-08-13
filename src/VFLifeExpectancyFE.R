# Clear workspace and load dependencies
setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
rm(list = ls())
DEBUG = FALSE
source("src/libs/loadLib.R")
library(dplyr)
library(sm)

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

dataIni = data %>% select(c(GDP = GDP.t0,LE = LE.t0, Year = Year.ini,countryCode))
dataFin = data %>% select(c(GDP = GDP.t1,LE = LE.t1, Year = Year.end,countryCode)) %>% filter(Year == max(Year))



data = rbind(dataIni,dataFin) %>% arrange(countryCode)

data$GDP = log(data$GDP)

nObs = length(unique(data$countryCode))
nT = length(unique(data$Year))
nEval = 2500


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

VF_hat <- panel_vf_results$estimator
X0_raw <- panel_vf_results$X0_raw

effects <- get_effects(panel_vf_results, X_obs = X0_raw, FE = TRUE, TE = TRUE)
FE_hat <- effects$alpha_i
TE_hat <- effects$gamma_t

bootstrap_samples <- bootstrapPanelVF(panel_vf_results, B = 100)

VF_bootstrap <- bootstrap_samples$estimators_array
FE_bootstrap <- bootstrap_samples$FE_array
TE_bootstrap <- bootstrap_samples$TE_array

VF_signif <- significanceBootstrap(VF_hat, VF_bootstrap, p_crit = 0.01)
FE_signif <- significanceBootstrap(FE_hat, FE_bootstrap, p_crit = 0.01)
TE_signif <- significanceBootstrap(TE_hat, TE_bootstrap, p_crit = 0.01)

xGrid = seq(from=min(X[,1,]), to=max(X[,1,]), length.out=round(sqrt(nEval)))
yGrid = seq(from=min(X[,2,]), to=max(X[,2,]), length.out=round(sqrt(nEval)))
x = as.matrix(expand.grid(xGrid, yGrid))

# Filter for significance
VF_hat[!VF_signif, ] <- 0

x[,1] = x[,1]/median(filter(data,Year==2010)$GDP)
x[,2] = x[,2]/median(filter(data,Year==2010)$LE)
VF_hat[,1] = VF_hat[,1]/median(filter(data,Year==2010)$GDP)
VF_hat[,2] = VF_hat[,2]/median(filter(data,Year==2010)$LE)


lengthArrows=(5/timeInterval)*1e-1*0.5
X0_raw_relative = X0_raw
X0_raw_relative[,1] = X0_raw_relative[,1]/median(filter(data,Year==2010)$GDP)
X0_raw_relative[,2] = X0_raw_relative[,2]/median(filter(data,Year==2010)$LE)
dev.new()
op <- par(family = "mono") #Possible families: "mono", "Helvetica","Palatino" or "Times" 
plot(x, type = "n", xlab = "Log GDP per capita (PPP in million 2017 USD, relative to 2010)", ylab = "Life Expectancy (relative to 2010)", main = "")
abline(h=1,lty=3)
abline(v=1,lty=3)
arrows(x[,1],x[,2],x[,1]+lengthArrows*VF_hat[,1],x[,2]+lengthArrows*VF_hat[,2],angle=15,col="black",length=0.05)
dataFirst = filter(data,Year == min(Year))
dataLatest = filter(data,Year == max(Year))
points(X0_raw_relative[,1],X0_raw_relative[,2],pch=19,col="red",cex=0.5)
est.dens = sm.density(X0_raw_relative,display="none")
contour(est.dens$eval.points[,1], est.dens$eval.points[,2], est.dens$estimate,add=T,col="purple")
dev.copy2pdf(file="outpics/VF0.pdf",width=7,height=7,family = "mono")

## plot FE ----
GDP0 = arrange(filter(data,Year == min(Year)),countryCode)$GDP
LE0 = arrange(filter(data,Year == min(Year)),countryCode)$LE
names0 = arrange(filter(data,Year == min(Year)),countryCode)$countryCode

cex_vals <- ifelse(FE_signif, 1.0, 0.1)

dev.new()
plot(GDP0,FE_hat[,1], main = "", xlab = "GDP (log)", ylab = "FE", pch=19, cex = cex_vals)
abline(h=0)
grid()
text(GDP0,FE_hat[,1],labels=names0,cex=0.5,pos=4,col="red")
dev.copy2pdf(file="outpics/FE_GDP0.pdf",width=7,height=7,family = "mono")

dev.new()
plot(LE0,FE_hat[,2], main = "", xlab = "LE", ylab = "FE", pch=19, cex = cex_vals)
abline(h=0)
grid()
text(LE0,FE_hat[,2],labels=names0,cex=0.5,pos=4,col="red")
dev.copy2pdf(file="outpics/FE_LE0.pdf",width=7,height=7,family = "mono")

## plot TE ----
times = sort(unique(data$Year))
cex_vals <- ifelse(TE_signif, 1.0, 0.1)

TE_quantiles <- apply(TE_bootstrap, c(1, 2), quantile, probs = c(0.025, 0.975), na.rm = TRUE)
lower_bound <- TE_quantiles[1, , ]
upper_bound <- TE_quantiles[2, , ]

dev.new()
plot(times[2:nT],TE_hat[,1], main = "", xlab = "Years", ylab = "TE GDP (log)", type="n", ylim = range(c(lower_bound[,1], upper_bound[,1], TE_hat[,1]), na.rm = TRUE))
polygon(c(times[2:nT], rev(times[2:nT])), c(lower_bound[, 1], rev(upper_bound[, 1])), col = "grey80", border = NA)
grid()
abline(h=0)
lines(times[2:nT],TE_hat[,1], main = "", xlab = "Years", ylab = "TE GDP (log)", pch=19, cex = cex_vals,type="b")
dev.copy2pdf(file="outpics/TE_GDP0.pdf",width=7,height=7,family = "mono")

dev.new()
plot(times[2:nT],TE_hat[,2], main = "", xlab = "Years", ylab = "TE LE", type="n", ylim = range(c(lower_bound[,2], upper_bound[,2], TE_hat[,2]), na.rm = TRUE))
polygon(c(times[2:nT], rev(times[2:nT])), c(lower_bound[, 2], rev(upper_bound[, 2])), col = "grey80", border = NA)
grid()
abline(h=0)
lines(times[2:nT],TE_hat[,2], main = "", xlab = "Years", ylab = "TE LE", pch=19, cex = cex_vals,type="b")
grid()
dev.copy2pdf(file="outpics/TE_LE0.pdf",width=7,height=7,family = "mono")


