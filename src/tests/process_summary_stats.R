# ==============================================================================
# SUMMARY STATISTICS PROCESSOR FOR PANEL DATA ESTIMATION TESTS
# ==============================================================================
#
# Purpose:
# This script processes summary statistics files generated from panel data
# estimation Monte Carlo simulations. It extracts performance metrics from text
# files and consolidates them into structured CSV files for further analysis.
#
# What it does:
# 1. Scans multiple test directories (e.g., test_pics_DistRot, test_pics_Rotation,
#    test_pics_DoubleWell) for summary statistics files
# 
# 2. Parses each summary text file to extract:
#    - Kernel type (Epanechnikov or Gaussian)
#    - Ridge regression parameter
#    - Condition number threshold
#    - Performance statistics for three components:
#      * Vector Field (VF) norm differences: mean, SD, RMSE, min, max
#      * Fixed Effects (FE) norm differences: mean, SD, RMSE, min, max
#      * Time Effects (TE) norm differences: mean, SD, RMSE, min, max
#
# 3. Organizes the extracted data by:
#    - Data Generating Process (DGP) type (extracted from directory name)
#    - Kernel type (Epa vs Gauss)
#
# 4. Outputs separate CSV files for each kernel type in each directory:
#    - Epa_summary.csv: Contains all results using Epanechnikov kernel
#    - Gauss_summary.csv: Contains all results using Gaussian kernel
#
# Input:
# - Text files matching pattern "Summary_Stats_Synthetic_*.txt" in specified
#   test directories, containing formatted statistics from simulation runs
#
# Output:
# - CSV files (Epa_summary.csv, Gauss_summary.csv) in each processed directory
#   with columns for all extracted parameters and statistics
#
# ==============================================================================

# Load necessary libraries
library(stringr)
library(dplyr)

# List of directories to process
data_dirs <- c("test_pics_DistRot", "test_pics_Rotation", "test_pics_DoubleWell")

# Function to parse a single summary file
parse_summary_file <- function(file_path, dgp_name) {
  # Extract parameters from the filename
  filename <- basename(file_path)
  
  kernel_match <- str_match(filename, "Synthetic_(Epa|Gauss)")
  kernel <- if (!is.na(kernel_match[1, 2])) kernel_match[1, 2] else NA
  
  ridge_match <- str_match(filename, "ridge_([\\de\\.\\-]+)_")
  ridge <- if (!is.na(ridge_match[1, 2])) as.numeric(ridge_match[1, 2]) else NA
  
  cond_match <- str_match(filename, "cond([\\de\\.\\-]+)\\.txt")
  condition_threshold <- if (!is.na(cond_match[1, 2])) as.numeric(cond_match[1, 2]) else NA
  
  # Read the file content
  lines <- readLines(file_path)
  
  # Function to extract stats from a line
  extract_stats <- function(line) {
    values <- as.numeric(str_extract_all(line, "[0-9\\.\\-]+")[[1]])
    if (length(values) >= 3) {
        # Mean, SD, RMSE
        return(list(mean = values[1], sd = values[2], rmse = values[3]))
    }
    return(list(mean = NA, sd = NA, rmse = NA))
  }
  
  extract_min_max <- function(line) {
    values <- as.numeric(str_extract_all(line, "[0-9\\.\\-]+")[[1]])
    if (length(values) == 2) {
        # Min, Max
        return(list(min = values[1], max = values[2]))
    }
     return(list(min = NA, max = NA))
  }

  # Find the lines with the statistics
  vf_stats_line <- lines[which(str_detect(lines, "VECTOR FIELD NORM DIFFERENCES")) + 2]
  vf_min_max_line <- lines[which(str_detect(lines, "VECTOR FIELD NORM DIFFERENCES")) + 3]
  fe_stats_line <- lines[which(str_detect(lines, "FIXED EFFECTS NORM DIFFERENCES")) + 2]
  fe_min_max_line <- lines[which(str_detect(lines, "FIXED EFFECTS NORM DIFFERENCES")) + 3]
  te_stats_line <- lines[which(str_detect(lines, "TIME EFFECTS NORM DIFFERENCES")) + 2]
  te_min_max_line <- lines[which(str_detect(lines, "TIME EFFECTS NORM DIFFERENCES")) + 3]
  
  # Extract the stats
  vf_stats <- extract_stats(vf_stats_line)
  vf_min_max <- extract_min_max(vf_min_max_line)
  fe_stats <- extract_stats(fe_stats_line)
  fe_min_max <- extract_min_max(fe_min_max_line)
  te_stats <- extract_stats(te_stats_line)
  te_min_max <- extract_min_max(te_min_max_line)
  
  # Combine into a single data frame row
  data_row <- data.frame(
    dgp = dgp_name,
    kernel = kernel,
    ridge = ridge,
    condition_threshold = condition_threshold,
    vf_mean = vf_stats$mean,
    vf_sd = vf_stats$sd,
    vf_rmse = vf_stats$rmse,
    vf_min = vf_min_max$min,
    vf_max = vf_min_max$max,
    fe_mean = fe_stats$mean,
    fe_sd = fe_stats$sd,
    fe_rmse = fe_stats$rmse,
    fe_min = fe_min_max$min,
    fe_max = fe_min_max$max,
    te_mean = te_stats$mean,
    te_sd = te_stats$sd,
    te_rmse = te_stats$rmse,
    te_min = te_min_max$min,
    te_max = te_min_max$max,
    source_file = filename
  )
  
  return(data_row)
}

# Loop over each directory and process the files
for (dir_name in data_dirs) {
  current_dir <- file.path(dir_name)
  
  # Get a list of all summary text files in the current directory
  files <- list.files(current_dir, pattern = "^Summary_Stats_Synthetic_.*\\.txt$", full.names = TRUE)
  
  if (length(files) == 0) {
    print(paste("No summary files found in", current_dir))
    next
  }
  
  # Extract the Data Generating Process name from the directory name
  dgp_name <- gsub("test_pics_", "", dir_name)

  # Initialize lists to store data for each kernel type
  epa_data <- list()
  gauss_data <- list()
  
  # Process each file and sort into the respective lists
  for (file in files) {
    parsed_data <- parse_summary_file(file, dgp_name)
    if (!is.na(parsed_data$kernel)) {
      if (parsed_data$kernel == "Epa") {
        epa_data[[length(epa_data) + 1]] <- parsed_data
      } else if (parsed_data$kernel == "Gauss") {
        gauss_data[[length(gauss_data) + 1]] <- parsed_data
      }
    }
  }
  
  # Combine lists into data frames and write to CSV
  if (length(epa_data) > 0) {
    epa_df <- bind_rows(epa_data)
    output_path <- file.path(current_dir, "Epa_summary.csv")
    write.csv(epa_df, output_path, row.names = FALSE)
    print(paste("Successfully created:", output_path))
  } else {
    print(paste("No Epa summary files found in", current_dir))
  }
  
  if (length(gauss_data) > 0) {
    gauss_df <- bind_rows(gauss_data)
    output_path <- file.path(current_dir, "Gauss_summary.csv")
    write.csv(gauss_df, output_path, row.names = FALSE)
    print(paste("Successfully created:", output_path))
  } else {
    print(paste("No Gauss summary files found in", current_dir))
  }
}
