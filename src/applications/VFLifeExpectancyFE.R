# Clear workspace and load dependencies
# setwd("~/Library/CloudStorage/OneDrive-UniversityofPisa/timeSpaceEvolutionEcAct/RVF/R code/Vector Field Estimation/")
rm(list = ls())
DEBUG = FALSE
source("src/libs/loadLib.R")
library(dplyr)
library(sm)
library(plm)

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
nTnEval = 2500

# Create panel data structure using plm
panel_data <- pdata.frame(data, index = c("countryCode", "Year"))



## qui tutto stima ----
# create array 
X = array(NA,dim=c(nObs, 2, nT))
for (year in unique(data$Year)) {
  X[,,(year - min(data$Year))/timeInterval+1] = cbind(data$GDP[data$Year == year], data$LE[data$Year == year])
}


kernels_to_test <- c("gauss", "epa")
ridge_params_to_test <- c(1e-2)
rcond_thresholds_to_test <- c(1e-3, 1e-4, 1e-6, 1e-8)

total_iterations <- length(kernels_to_test) * length(ridge_params_to_test) * length(rcond_thresholds_to_test)
current_iteration <- 0

for (kernel in kernels_to_test) {
  for (ridge in ridge_params_to_test) {
    for (rcond in rcond_thresholds_to_test) {
      current_iteration <- current_iteration + 1
      
      kernel_type_str <- tools::toTitleCase(kernel)
      ridge_str <- if (ridge == 0) "0" else gsub("e-0", "e-", formatC(ridge, format = "e", digits = 0))
      rcond_str <- if (rcond == 0) "0" else gsub("e-0", "e-", formatC(rcond, format = "e", digits = 0))
      filename <- sprintf("%s_ridge_%s_cond%s.svg", kernel_type_str, ridge_str, rcond_str)
      
      cat(sprintf("\n--- [ %d / %d ] Processing: %s ---\n", current_iteration, total_iterations, filename))

      analysis_results <- runPanelVFAnalysis(panel_data,
                                             var_cols = c("GDP", "LE"),
                                             id_col = "countryCode",
                                             time_col = "Year",
                                             nEval = nEval,
                                             FE = TRUE,
                                             TE = TRUE,
                                             uniform_weights = TRUE,
                                             kernel.type = kernel,
                                             method.h = "silverman",
                                             chunk_size = 1024,
                                             bootstrap_B = 100, 
                                             ridge_param = ridge,
                                             rcond_threshold = rcond)
      
      plotPanelVFAnalysis(analysis_results,
                          timeInterval = timeInterval,
                          rescale = TRUE,
                          rescale_ref_index = nT,
                          years = sort(unique(data$Year)),
                          label_names = arrange(filter(data,Year == min(Year)),countryCode)$countryCode,
                          component_names = c("Log GDP", "Life Expectancy"),
                          save_plots = TRUE)
    }
  }
}
