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

# Call the new encapsulated function from panel.R
analysis_results <- runPanelVFAnalysis(X,
                                     nEval = nEval,
                                     FE = TRUE,
                                     TE = TRUE,
                                     uniform_weights = TRUE,
                                     kernel.type = "epa",
                                     method.h = "silverman",
                                     chunk_size = 512,
                                     bootstrap_B = 100)

plotPanelVFAnalysis(analysis_results,
                    timeInterval = timeInterval,
                    rescale = TRUE,
                    rescale_ref_index = nT,
                    years = sort(unique(data$Year)),
                    label_names = arrange(filter(data,Year == min(Year)),countryCode)$countryCode,
                    component_names = c("Log GDP", "Life Expectancy"),
                    out_dir = "outpics",
                    save_plots = TRUE,
                    show_plots = TRUE)

