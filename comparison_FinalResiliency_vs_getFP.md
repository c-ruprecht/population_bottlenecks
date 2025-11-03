# Comparison: FinalResiliencyWithCFU.R (getNrNb) vs getFP_12112023.R (getFP)

## Overview
- **FinalResiliencyWithCFU.R**: Implements the "resiliency" method as described in the original paper - iteratively removes barcodes and finds plateau values
- **getFP_12112023.R**: Uses a simulation-based approach for founding population estimation

---

## Key Output Parameters

| Parameter | FinalResiliencyWithCFU.R | getFP_12112023.R |
|-----------|-------------------------|------------------|
| **Primary outputs** | Nb, Nr | Nb, Ns, Ns_MinCutoff |
| **Terminology** | Nr = "resilient" founding population | Ns = founding population size |
| **Column names** | "Nb", "Nr" (+ calibrated versions if provided) | "TotalReads", "Number of barcodes", "Ns_MinCutoff", "Nb", "Ns", "AverageFrequency", "CFU", "Log10Ns", "Log10CFU", "CFU/Ns" |

---

## Detailed Parameter Comparison

### 1. **Nb (Bottleneck Size)**

#### FinalResiliencyWithCFU.R (Lines 68-82)
```r
# Iteratively calculated - removes barcodes one at a time
minusone <- function() {
  # Wright's F-statistic
  F <- mean(na.omit(sigma))
  nb <- 1 / (F - 1/sum(invec) - 1/sum(outvec))
  x <<- c(x, nb)  # Append to vector
  bindsorted <<- bindsorted[1:dim(bindsorted)[1]-1,]  # Remove most abundant barcode
}
replicate(times, minusone())  # Creates resiliency curve
```
- **Method**: Iteratively removes most abundant barcodes, calculates Nb after each removal
- **Result**: Vector `x` contains Nb at each iteration (the resiliency curve)
- **Final Nb**: `x[1]` - the first value (Nb before any removals)

#### getFP_12112023.R (Lines 105-116, 243-250)
```r
# Calculated once for subsets, not iteratively
minusone <- function(dims) {
  input <- bindsorted[,1][1:dims]  # Take subset
  out <- bindsorted[,2][1:dims]
  # ... Wright's F-statistic on subset
  max(1 / (F - 1/sum(invec) - 1/sum(outvec)), 1/F)
}
x <- sapply(timestorep, minusone)  # Calculate for different subset sizes

# Also calculated once on full cleaned dataset (lines 243-250)
Nb <- 1 / (F - 1/sum(invec) - 1/sum(outvec))
```
- **Method**: Calculates Nb for pre-defined subset sizes (not iterative removal)
- **Result**: Vector `x` contains Nb for different numbers of included barcodes
- **Final Nb**: Single value calculated on full cleaned dataset after noise removal

**Key Difference**: FinalResiliencyWithCFU iteratively removes; getFP calculates for subsets

---

### 2. **Nr/Ns (Founding Population Size)**

#### FinalResiliencyWithCFU.R - Nr (Lines 210-239)
```r
# After noise removal, re-run resiliency on cleaned data
z <- 0
minusonefinal <- function() {
  # Same Wright's F-statistic on noise-filtered data
  nb <- 1 / (F - 1/sum(invec) - 1/sum(outvec))
  z <<- c(z, nb)
  bindsortedcopy2 <<- bindsortedcopy2[1:dim(bindsortedcopy2)[1]-1,]
}
replicate(sum(bindsortedcopy2[2] != 0), minusonefinal())

# Nr is the MAXIMUM Nb from the cleaned resiliency curve
finalvalue <- max(z)
if(finalvalue < x[1]) {finalvalue <- x[1]}  # Ensure Nr >= Nb

# Additional constraint from barcode count simulation
minimumnearesttimes <- max(indices[,3][(indices[,3] <= noisestart)])
CompBottleneck <- NLSstClosestX(dfxy, minimumnearesttimes)
if(finalvalue < CompBottleneck) {finalvalue <- CompBottleneck}
```
- **Method**: Runs resiliency again after noise filtering, takes max plateau value
- **Basis**: Maximum Nb value from iterative resiliency curve on cleaned data
- **Constraint**: Must be >= original Nb and >= CompBottleneck (from barcode simulation)
- **Philosophy**: The "resilient" bottleneck estimate (plateau after removing abundant clones)

#### getFP_12112023.R - Ns (Lines 221-241)
```r
# Resampling simulation approach
FirstResample <- rmultinom(1, sum(outvecwithoutnoise), ReferenceVector/sum(ReferenceVector))

GetNewBotTable <- function(n) {
  vec <- rmultinom(5, n, FirstResample/sum(FirstResample))
  vec[vec!=0] <- 1
  mean(colSums(vec))  # How many barcodes detected from n cells?
}

steps1 <- round(seq(from = 1, to = sum(invec != 0), length.out = 100))
steps2 <- round(seq(from = sum(invec != 0), to = sum(invec != 0)*20, length.out = 100))
steps <- c(steps1, steps2)
y <- as.numeric(sapply(steps, GetNewBotTable))

# Create simulation curve: x=population size, y=detected barcodes
dfxy2 <- floor(sortedXyData(x = xvals, y = yvals))

# Find population size that would produce observed barcode count
Ns <- NLSstClosestX(dfxy2, noisestart)
```
- **Method**: Simulates how many barcodes would be detected from different founding population sizes
- **Basis**: Creates theoretical curve, finds x-value (population size) matching observed barcode count
- **Philosophy**: Statistical inference from barcode diversity via resampling theory
- **Does NOT use resiliency plateau method**

**Key Difference**:
- **Nr**: Max value from iterative Nb calculations (empirical plateau)
- **Ns**: Population size estimate from simulation curve (theoretical inference)

---

### 3. **Ns_MinCutoff / Ns_MinWeight**

#### FinalResiliencyWithCFU.R - Not directly calculated
- The script uses `minweight` for noise detection logic (lines 180-183)
- Does NOT output a separate "Ns at minweight" metric
- `minweight` affects where noise starts but no explicit minweight-based founding population estimate

#### getFP_12112023.R - Ns_MinCutoff (Lines 254-257)
```r
# Calculate which barcodes account for (1-minweight) of reads
ReadsAtWeightCutoff <- (1-minweight)*sum(outvec)
MinCutoff <- outvecsorted[min(which(cumsum(outvecsorted[outvecsorted > 0]) > ReadsAtWeightCutoff))]

# Count barcodes above this frequency threshold
NBarcodesAtMinweight <- sum(outvec > MinCutoff) + 1

# Use simulation curve to estimate Ns for this barcode count
Ns_MinCutoff <- NLSstClosestX(dfxy2, NBarcodesAtMinweight)
```
- **Method**: Threshold-based approach using cumulative frequency
- **Basis**: Counts barcodes that collectively account for ≥(1-minweight)% of reads
- **Philosophy**: Conservative estimate less dependent on algorithmic noise detection

**Key Difference**: Only getFP outputs this metric; FinalResiliencyWithCFU doesn't calculate it

---

### 4. **Noise Detection Method**

#### Both Scripts (Similar Logic)
Both use log-weight differences to find noise threshold:

**FinalResiliencyWithCFU.R (Lines 170-183)**
```r
weights <- log(indices[,2])  # Log of fraction of reads
weightsforsubtraction <- c(weights[2:length(values)], 0)
weightsdif <- (weights - weightsforsubtraction)[1:length(weights)-1]
noisestart <- indices[,3][which(weightsdif == max(weightsdif))]

# Adjust if no populations below minweight
if (is.na(indices[,3][min(which(indices[,2] < minweight))])) {
  noisestart <- max(indices[,3])
}
# Adjust if populations after noise sum to > minweight
if (sum(indices[,2][which(indices[,3] > noisestart)]) > minweight) {
  noisestart <- indices[,3][locationofminweightcutoff]
}
```

**getFP_12112023.R (Lines 194-206)**
```r
weights <- log(indices[,2])
weightsforsubtraction <- c(weights[2:length(values)], 0)
weightsdif <- (weights - weightsforsubtraction)[1:length(weights)-1]
noisestart <- indices[,1][which(weightsdif == max(weightsdif))]

# if no populations below minweight threshold
if (min(indices[,2]) > minweight) {noisestart <- max(indices[,1])}
# if populations after noise > minweight
if (sum(indices[,2][which(indices[,1] > noisestart)]) > minweight) {
  noisestart <- indices[,1][locationofminweightcutoff]
}
```

**Similarity**: Nearly identical noise detection logic using log-weight differences

---

### 5. **Population Structure Detection (Indices Table)**

#### FinalResiliencyWithCFU.R (Lines 147-168)
```r
# After finding local minima in resiliency curve
accountsfor <- function(t) {
  outvecsorted <- sort(outvec)
  topnumbers <- tail(outvecsorted, t)
  sum(topnumbers) / sum(outvec)
}

# Get max Nb up to each break point
getintervalnb <- function(y) {
  max(x[1:y])  # Max Nb in this segment
}
nbintervals <- unlist(lapply(guessesuniquesorted, getintervalnb))

indices <- data.frame(nbintervals, subtracted, guessesuniquesorted)
colnames(indices) <- c("Nr", "Accounts for", "Number of barcodes")
```
- **Output**: Each population segment gets an Nr value (max Nb in that segment)
- **Column name**: "Nr" (representing max bottleneck for that population)

#### getFP_12112023.R (Lines 177-191)
```r
accountsfor <- function(t) {
  outvecsorted <- sort(outvec)
  topnumbers <- tail(outvecsorted, t)
  sum(topnumbers) / sum(outvec)
}

fractionaccounted <- sapply(guessesuniquesorted, accountsfor)
staggered <- c(0, fractionaccounted)[1:length(fractionaccounted)]
subtracted <- fractionaccounted - staggered

indices <- data.frame(guessesuniquesorted, subtracted)
colnames(indices) <- c("Number of barcodes", "Fraction of reads")
```
- **Output**: Simple table of barcode counts and read fractions
- **Does NOT include Nb/Nr values for each segment**

**Key Difference**: FinalResiliencyWithCFU includes Nr for each population segment; getFP doesn't

---

### 6. **Calibration Support**

#### FinalResiliencyWithCFU.R (Lines 252-262, 278-279)
```r
if(!is.null(CalibrationFile)) {
  LookUpTable <- read.csv(CalibrationFile, row.names = 1)
  nbcalibrated <- as.numeric(LookUpTable[as.character(round(log10(x[1]), digits = 2)),])
  nrcalibrated <- as.numeric(LookUpTable[as.character(round(log10(finalvalue), digits = 2)),])
}

# Output columns
if(!is.null(CalibrationFile)) {
  colnames(TableOfEstimates) <- c("Nb", "Nr", "NbLower_Cal", "NbMedian_Cal",
                                  "NbUpper_Cal", "NrLower_Cal", "NrMedian_Cal", "NrUpper_Cal")
}
```
- **Supports**: Optional calibration lookup table
- **Output**: 3 calibrated values each for Nb and Nr (Lower, Median, Upper)

#### getFP_12112023.R
- **No calibration support**

---

### 7. **Additional Outputs**

#### FinalResiliencyWithCFU.R Output Columns
```
Without calibration: "Nb", "Nr", "Nb_cal", "Nr_cal"
With calibration:    "Nb", "Nr", "NbLower_Cal", "NbMedian_Cal", "NbUpper_Cal",
                     "NrLower_Cal", "NrMedian_Cal", "NrUpper_Cal"
```

#### getFP_12112023.R Output Columns
```
"TotalReads", "Number of barcodes", "Ns_MinCutoff", "Nb", "Ns",
"AverageFrequency", "CFU", "Log10Ns", "Log10CFU", "CFU/Ns"
```

**Additional metrics in getFP only:**
- **TotalReads**: Sum of all reads in sample
- **AverageFrequency**: `1/geoMean(frequencies)` - geometric mean-based diversity metric
- **Log10Ns, Log10CFU**: Log-transformed values
- **CFU/Ns**: Ratio of CFU to founding population

---

### 8. **InocCFU Parameter**

#### FinalResiliencyWithCFU.R (Lines 18-29)
```r
# Uses InocCFU to create theoretical barcode detection curve
RoundedRefVector <- round(round(ReferenceVector) * InocCFU / sum(round(ReferenceVector)))

GetBotTable <- function(n) {
  vec <- as.numeric(rmvhyper(1, RoundedRefVector, n))  # Multivariate hypergeometric
  sum(vec!=0)
}
```
- **Purpose**: Creates simulation curve for estimating Nr from barcode counts
- **Method**: Uses multivariate hypergeometric distribution (`rmvhyper` from `extraDistr`)
- **Required**: Yes, must specify inoculum CFU

#### getFP_12112023.R (Lines 221-233)
```r
# Does NOT use InocCFU parameter
FirstResample <- rmultinom(1, sum(outvecwithoutnoise), ReferenceVector/sum(ReferenceVector))

GetNewBotTable <- function(n) {
  vec <- rmultinom(5, n, FirstResample/sum(FirstResample))  # Multinomial
  vec[vec!=0] <- 1
  mean(colSums(vec))
}
```
- **Purpose**: Creates simulation curve using only observed data
- **Method**: Uses multinomial distribution (`rmultinom`)
- **Required**: No InocCFU parameter

**Key Difference**: FinalResiliencyWithCFU requires experimental InocCFU; getFP doesn't

---

### 9. **Noise Correction**

#### FinalResiliencyWithCFU.R (Lines 42-47)
```r
if (CorrectForNoise != 0) {
  ResampledInvec <- rmvhyper(1, round(invec), round(CorrectForNoise*sum(outvec)))
  outvec <- as.numeric(outvec - ResampledInvec)
  outvec[outvec < 0] <- 0
}
```
- **Parameter**: `CorrectForNoise` (numeric, e.g., 0.001)
- **Method**: Subtracts expected noise using hypergeometric resampling
- **Applied**: Once at the beginning, before any analysis

#### getFP_12112023.R (Lines 48-59)
```r
# Complex noise adjustment based on singleton ratio
oneratio <- sqrt(sum(outvec > 1) / sum(outvec == 1))
if(oneratio < .8) {oneratio <- -(oneratio^2)}
NoiseCorrectFactor <- 10^-oneratio
NoiseDist <- rowSums(ReadsTable[,-WhereAreReferences])

ResampledInvec <- rmultinom(50, round(sum(outvec) * NoiseCorrectFactor), NoiseDist)
ResampledInvec[ResampledInvec > .005*sum(outvec)] <- 0
outvec <- ceiling(rowMeans(outvec - ResampledInvec))
outvec[outvec < 0] <- 0
```
- **Parameter**: None (automatic)
- **Method**: Calculates noise correction factor from singleton/non-singleton ratio
- **Applied**: Automatically using adaptive noise model

**Key Difference**: FinalResiliencyWithCFU uses user-specified correction; getFP uses automatic

---

## Summary Table

| Feature | FinalResiliencyWithCFU.R | getFP_12112023.R |
|---------|-------------------------|------------------|
| **Core method** | Iterative resiliency (paper method) | Simulation-based inference |
| **Nb calculation** | Iterative removal, x[1] | Single calculation on full data |
| **Founding pop metric** | Nr (max from resiliency curve) | Ns (from simulation curve) |
| **Ns_MinCutoff** | ❌ Not calculated | ✅ Calculated |
| **Requires InocCFU** | ✅ Yes (required parameter) | ❌ No |
| **Noise correction** | User-specified (CorrectForNoise) | Automatic (singleton-based) |
| **Calibration support** | ✅ Yes (optional) | ❌ No |
| **Population segments** | Includes Nr for each | Only counts/fractions |
| **Additional metrics** | Minimal | TotalReads, AverageFreq, CFU ratios |
| **Dependencies** | `extraDistr`, `EnvStats` | `EnvStats` only |
| **Matches paper** | ✅ Yes | ❌ No (different algorithm) |

---

## Which to Use?

**Use FinalResiliencyWithCFU.R (getNrNb) when:**
- You want the method described in the published paper
- You have known inoculum CFU
- You want calibrated estimates
- You need Nr (resilient bottleneck estimate)

**Use getFP_12112023.R (getFP) when:**
- You don't know inoculum CFU
- You want Ns_MinCutoff (frequency threshold-based estimate)
- You prefer automatic noise correction
- You want additional summary statistics (AverageFrequency, CFU/Ns ratios)
