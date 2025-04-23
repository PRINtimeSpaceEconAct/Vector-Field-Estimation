#' Loads all required libraries and source files for vector field estimation
#' 
#' This script loads the necessary R packages and sources all the local
#' library files needed for vector field estimation. It should be called
#' before using any of the vector field estimation functions.
#'
library(extremefit)
library(Matrix)
library(spatstat)
library(provenance)
library(MASS)
library(expm)
library(emdbook)
library(pracma)

source("src/libs/utils.R")
source("src/libs/prettyprinting.R")
source("src/libs/distances.R")
source("src/libs/kernel.R")
source("src/libs/density.R")
source("src/libs/NW.R")
source("src/libs/localLinear.R")
source("src/libs/Variance.R")
source("src/libs/Forecast.R")
