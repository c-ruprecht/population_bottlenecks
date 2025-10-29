#!/usr/bin/env Rscript

# Load the getFP function
source("test/holmes2024/getFP_12112023.R")

# Run getFP on subset
# Reference columns are 1-6 in R (1-indexed)
# Sample is s43_ch06 (column 7)

getFP(
  ReadsTableName = "test/holmes2024/inputs/NoHopFreq_Masked_tatC_subset.csv",
  CFUtable = "test/holmes2024/inputs/CFU_tatC_12112023.csv",
  WhereAreReferences = c(1, 2, 3, 4, 5, 6),
  minweight = 0.03,
  outputfilename = "test/holmes2024/r_subset_output/NsNb_subset.csv"
)

# Print results
print("R getFP completed successfully")
print(TableOfEstimates)
