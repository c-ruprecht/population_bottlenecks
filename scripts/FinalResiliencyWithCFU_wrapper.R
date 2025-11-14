#!/usr/bin/env Rscript

# Wrapper script for FinalResiliencyWithCFU.R to enable command-line usage
# Usage: Rscript FinalResiliencyWithCFU_wrapper.R <ReadsTable> <CFUtable> <InocCFU> <WhereAreReferences> <minweight> <CorrectForNoise> <output_dir> [CalibrationFile]

# Load required packages
# Note: These packages must be pre-installed in your user library
# To install them, run on the login node:
#   ml R/4.2.0
#   R
#   install.packages(c("extraDistr", "EnvStats"))
library(extraDistr)
library(EnvStats)

# Source the error handling version of getNrNb function
source("/hpc/users/ruprec02/git/population_bottlenecks/scripts/FinalResiliencyWithCFU_ErrorHandling.R")

# Parse command-line arguments
cmd_args <- commandArgs(trailingOnly = TRUE)

if (length(cmd_args) < 7 || length(cmd_args) > 8) {
  stop("Usage: Rscript FinalResiliencyWithCFU_wrapper.R <ReadsTable> <CFUtable> <InocCFU> <WhereAreReferences> <minweight> <CorrectForNoise> <output_dir> [CalibrationFile]")
}

ReadsTable <- cmd_args[1]
CFUtable <- cmd_args[2]
if (CFUtable == "NULL") {
  CFUtable <- NULL
}
InocCFU <- as.numeric(cmd_args[3])
WhereAreReferences <- as.integer(unlist(strsplit(cmd_args[4], ",")))
minweight <- as.numeric(cmd_args[5])
CorrectForNoise <- as.numeric(cmd_args[6])
output_dir <- cmd_args[7]

# Optional calibration file
CalibrationFile <- NULL
if (length(cmd_args) == 8 && cmd_args[8] != "NULL") {
  CalibrationFile <- cmd_args[8]
}

# Print parsed arguments for debugging
print("Parsed arguments:")
print(paste("ReadsTable:", ReadsTable))
print(paste("CFUtable:", CFUtable))
print(paste("InocCFU:", InocCFU))
print(paste("WhereAreReferences:", paste(WhereAreReferences, collapse=",")))
print(paste("minweight:", minweight))
print(paste("CorrectForNoise:", CorrectForNoise))
print(paste("output_basename:", output_dir))
print(paste("CalibrationFile:", ifelse(is.null(CalibrationFile), "NULL", CalibrationFile)))

# The output_dir parameter is actually a basename (e.g., /path/to/FP/experiment/strain)
# The R function will append _TableOfEstimates.csv and _FrequenciesWithoutNoise.csv
outputfilename <- output_dir  # This is the basename prefix

# Create the output directory
output_directory <- dirname(outputfilename)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# Call the getNrNb function with the basename
# The function will create:
#   outputfilename + '_TableOfEstimates.csv'
#   outputfilename + '_FrequenciesWithoutNoise.csv'
getNrNb(ReadsTable, CFUtable, InocCFU, WhereAreReferences, minweight, CorrectForNoise, outputfilename, CalibrationFile)

print(paste("Analysis complete. Output written to:", output_directory))
print("Files created:")
print(paste("  -", paste0(outputfilename, "_TableOfEstimates.csv")))
print(paste("  -", paste0(outputfilename, "_FrequenciesWithoutNoise.csv")))
