#!/usr/bin/env Rscript

# Wrapper script for getFP_12112023.R to enable command-line usage
# Usage: Rscript getFP_12112023_wrapper.R <ReadsTable> <CFUtable> <WhereAreReferences> <minweight> <output_dir>

# Source the original getFP function
source("/hpc/users/ruprec02/git/population_bottlenecks/test/holmes2024/getFP_12112023.R")

# Parse command-line arguments
cmd_args <- commandArgs(trailingOnly = TRUE)

if (length(cmd_args) != 5) {
  stop("Usage: Rscript getFP_12112023_wrapper.R <ReadsTable> <CFUtable> <WhereAreReferences> <minweight> <output_dir>")
}

ReadsTable <- cmd_args[1]
CFUtable <- cmd_args[2]
if (CFUtable == "NULL") {
  CFUtable <- NULL
}
WhereAreReferences <- as.integer(unlist(strsplit(cmd_args[3], ",")))
minweight <- as.numeric(cmd_args[4])
output_dir <- cmd_args[5]

# Print parsed arguments for debugging
print("Parsed arguments:")
print(paste("ReadsTable:", ReadsTable))
print(paste("CFUtable:", CFUtable))
print(paste("WhereAreReferences:", WhereAreReferences))
print(paste("minweight:", minweight))
print(paste("output_dir:", output_dir))

# Convert output_dir to outputfilename format expected by getFP_12112023.R
# The old function expects a single CSV output filename, not a directory
# We'll create the directory structure and pass the filename
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
outputfilename <- file.path(output_dir, "TableOfEstimates.csv")

# Call the getFP function with plots=FALSE
# Note: The original getFP_12112023.R doesn't expose plots parameter to the main function,
# but ResiliencyIndices has plots=FALSE as default, so we're good
getFP(ReadsTable, CFUtable, WhereAreReferences, minweight, outputfilename)

# The function writes:
# - TableOfEstimates.csv (using outputfilename parameter)
# - FrequenciesWithoutNoise.csv (using same directory as outputfilename)
# Both files will be in output_dir

print(paste("Analysis complete. Output written to:", output_dir))
