#!/usr/bin/env Rscript

# Load the getFP function
source("test/holmes2024/getFP_12112023.R")

# Run getFP on 3-sample subset
# Reference columns are 1-6 in R (1-indexed)
# Samples are s43_ch06, s44_ch06, s45_ch06 (columns 7-9)

getFP(
  ReadsTableName = "test/holmes2024/inputs/NoHopFreq_Masked_tatC_3samples.csv",
  CFUtable = "test/holmes2024/inputs/CFU_tatC_12112023.csv",
  WhereAreReferences = c(1, 2, 3, 4, 5, 6),
  minweight = 0.03,
  outputfilename = "test/holmes2024/r_3samples_output/NsNb_3samples.csv"
)

# Print results
print("R getFP completed successfully")
print("TableOfEstimates:")
print(TableOfEstimates)
