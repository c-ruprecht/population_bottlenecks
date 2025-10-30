#!/usr/bin/env Rscript

# Wrapper script for MajorityDistance_12112023.R to enable command-line usage
# Usage: Rscript MajorityDistance_12112023_wrapper.R <ReadsTableName> <expname> <ComparisonMetaDataName> <limit> <output_dir>

# Source the original ResiliencyGeneticDistance function
source("/hpc/users/ruprec02/git/population_bottlenecks/scripts/MajorityDistance_12112023.R")

# Parse command-line arguments
cmd_args <- commandArgs(trailingOnly = TRUE)

if (length(cmd_args) != 5) {
  stop("Usage: Rscript MajorityDistance_12112023_wrapper.R <ReadsTableName> <expname> <ComparisonMetaDataName> <limit> <output_dir>")
}

ReadsTableName <- cmd_args[1]
expname <- cmd_args[2]
ComparisonMetaDataName <- cmd_args[3]
limit <- as.numeric(cmd_args[4])
output_dir <- cmd_args[5]

# Print parsed arguments for debugging
print("Parsed arguments:")
print(paste("ReadsTableName:", ReadsTableName))
print(paste("expname:", expname))
print(paste("ComparisonMetaDataName:", ComparisonMetaDataName))
print(paste("limit:", limit))
print(paste("output_dir:", output_dir))

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Change working directory to output_dir so CSV files are written there
# The original function writes CSV files to the current working directory
original_wd <- getwd()
setwd(output_dir)

# Call the ResiliencyGeneticDistance function
tryCatch({
  ResiliencyGeneticDistance(ReadsTableName, expname, ComparisonMetaDataName, limit)
}, finally = {
  # Always restore the original working directory
  setwd(original_wd)
})

# The function writes:
# - CorrectedRD_<expname>.csv
# - GD_<expname>.csv
# - CorrectedFRD_<expname>.csv
# All files will be in output_dir

print(paste("Analysis complete. Output written to:", output_dir))
print(paste("Files created:"))
print(paste("  -", file.path(output_dir, paste0("CorrectedRD_", expname, ".csv"))))
print(paste("  -", file.path(output_dir, paste0("GD_", expname, ".csv"))))
print(paste("  -", file.path(output_dir, paste0("CorrectedFRD_", expname, ".csv"))))
