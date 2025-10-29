source("test/holmes2024/getFP_debug.R")

cat("\n\n========== Running debug getFP on s43_ch06 ==========\n\n")

result <- getFP(
  ReadsTableName = "test/holmes2024/inputs/NoHopFreq_Masked_tatC_3samples.csv",
  CFUtable = "test/holmes2024/inputs/CFU_tatC_12112023.csv",
  WhereAreReferences = c(1, 2, 3, 4, 5, 6),
  minweight = 0.03,
  outputfilename = "test/holmes2024/r_debug_output/NsNb_debug.csv"
)

cat("\n\n========== getFP complete ==========\n\n")
