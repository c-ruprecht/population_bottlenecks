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

# Source the original getNrNb function
source("/hpc/users/ruprec02/git/population_bottlenecks/scripts/FinalResiliencyWithCFU.R")

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
print(paste("output_dir:", output_dir))
print(paste("CalibrationFile:", ifelse(is.null(CalibrationFile), "NULL", CalibrationFile)))

# Convert output_dir to outputfilename format expected by FinalResiliencyWithCFU.R
# The function expects a single CSV output filename, not a directory
# We'll create the directory structure and pass the filename
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
outputfilename <- file.path(output_dir, "TableOfEstimates.csv")

# Change working directory to output_dir so FrequenciesWithoutNoise.csv is written there
# The original function writes FrequenciesWithoutNoise.csv to the current working directory
original_wd <- getwd()
setwd(output_dir)

# Call the getNrNb function
tryCatch({
  getNrNb(ReadsTable, CFUtable, InocCFU, WhereAreReferences, minweight, CorrectForNoise, outputfilename, CalibrationFile)
}, finally = {
  # Always restore the original working directory
  setwd(original_wd)
})

# The function writes:
# - TableOfEstimates.csv (using outputfilename parameter in output_dir)
# - FrequenciesWithoutNoise.csv (to current working directory, which we set to output_dir)
# Both files will be in output_dir

print(paste("Analysis complete. Output written to:", output_dir))
print("Files created:")
print(paste("  -", file.path(output_dir, "TableOfEstimates.csv")))
print(paste("  -", file.path(output_dir, "FrequenciesWithoutNoise.csv")))
