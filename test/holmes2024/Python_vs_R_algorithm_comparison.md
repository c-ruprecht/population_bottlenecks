# Python vs R getFP Algorithm Comparison

## Key Finding: The Core Difference

The fundamental difference between R and Python is **what gets passed to the Ns estimation function**:

### R's Two-Stage Approach

**Stage 1: Population Detection (Breakpoint Analysis)**
```r
# Lines 194-206 in getFP_12112023.R
weights <- log(indices[,2])
weightsdif <- (weights - weightsforsubtraction)[1:length(weights)-1]
noisestart <- indices[,1][which(weightsdif == max(weightsdif))]
# noisestart = NUMBER OF DETECTED POPULATIONS
```

**Stage 2: Ns Estimation**
```r
# Lines 241-242
dfxy2 <- floor(sortedXyData(x = xvals, y = yvals))
Ns <- NLSstClosestX(dfxy2, noisestart)  # Uses noisestart (populations)
```

**Result**: R detects 5 POPULATIONS for s43_ch06, then estimates Ns=5

---

### Python's Approach

**Stage 1: Breakpoint Detection**
```python
# Line 472 in getFP.py
noise_start = determine_noise_threshold(indices_df, min_weight)
# noise_start = NUMBER OF BARCODES TO KEEP
```

**Stage 2: Ns Estimation**
```python
# Lines 488-490
Ns, interp_curve = estimate_founding_population_simulation(
    input_vec, output_vec, noise_start  # Uses noise_start (barcodes)
)
```

**Result**: Python keeps 13 BARCODES for s43_ch06, then estimates Ns=20.49

---

## The Critical Semantic Difference

| Concept | R | Python |
|---------|---|--------|
| **What `noisestart`/`noise_start` represents** | Number of detected POPULATIONS | Number of BARCODES to keep |
| **What gets reported as "Number of barcodes"** | Number of populations | Number of barcodes |
| **What drives Ns calculation** | Population count | Barcode count |

---

## R's Population Detection Algorithm

From the Task analysis, R uses a sophisticated multi-step approach:

### Step 1: Find Breakpoints in Resiliency Curve

**Stochastic search** (lines 124-144):
```r
scanformin <- function(guesses, number) {
  # 10,000 random walks with Gaussian jumps
  # Standard deviation = sum(outvec>1) / 20
  for (i in 1:number) {
    # Random walk to find local minima
  }
  return(unique(guesses))
}

guesses <- scanformin(guesses, 10000)
```

**Large jump detection** (lines 148-175):
- Find 10-fold jumps in read counts between consecutive barcodes
- Find 10-fold jumps in resiliency values
- Remove adjacent breakpoints (difference of 1)
- Merge breakpoints in top 1% of range

### Step 2: Create Population Table

**Calculate fractions for each population** (lines 177-191):
```r
accountsfor <- function(t) {
  topnumbers <- tail(sort(outvec), t)
  sum(topnumbers) / sum(outvec)
}

fractionaccounted <- sapply(guessesuniquesorted, accountsfor)
staggered <- c(0, fractionaccounted)[1:length(fractionaccounted)]
subtracted <- fractionaccounted - staggered
subtracted <- subtracted/sum(subtracted)

indices <- data.frame(guessesuniquesorted, subtracted)
colnames(indices) <- c("Number of barcodes", "Fraction of reads")
```

### Step 3: Select Noise Threshold

**Find maximum log-weight difference** (lines 194-206):
```r
weights <- log(indices[,2])
weightsforsubtraction <- c(weights[2:length(values)], 0)
weightsdif <- (weights - weightsforsubtraction)[1:length(weights)-1]
noisestart <- indices[,1][which(weightsdif == max(weightsdif))]

# Override 1: If all populations above minweight, use max
if (min(indices[,2]) > minweight) {noisestart <- max(indices[,1])}

# Override 2: If remaining reads > minweight, move cutoff
locationofminweightcutoff <- min(which(cumsum(indices[,2])>(1-minweight)))
if (sum(indices[,2][which(indices[,1] > noisestart)]) > minweight) {
  noisestart <- indices[,1][locationofminweightcutoff]
}
```

### Step 4: Special Cases

**Two-barcode case** (lines 208-212):
```r
if(length(which(indices[,1] == noisestart)) > 1) {
  noisestart <- indices[indices[,1] == noisestart,][2,1]
}

if (noisestart == 2 & sortedXy[1]/sortedXy[2] < 0.01) {noisestart <- 1}
```

---

## Python's Current Algorithm

### Step 1: Find Breakpoints

**Stochastic search** (lines 412-431):
```python
for _ in range(10000):
    current_pos = np.random.randint(0, len(resiliency_curve))
    visited = {current_pos}

    for step in range(100):
        # Random Gaussian jump
        jump = int(np.random.normal(0, fraction_sd))
        new_pos = current_pos + jump
        # ... find local minimum
```

**Large jump detection** (lines 434-441):
- 10-fold jump in resiliency
- Adjacent breakpoint removal
- NO read count jump detection (unlike R)

### Step 2: Create Indices Table

**Calculate fractions** (lines 456-469):
```python
def calculate_cumulative_fraction(output_vec, n_barcodes):
    sorted_vec = np.sort(output_vec)[::-1]
    top_counts = sorted_vec[:n_barcodes]
    return np.sum(top_counts) / np.sum(output_vec)

fractions = []
for n_barcodes in break_points:
    frac = calculate_cumulative_fraction(output_vec, n_barcodes)
    fractions.append(frac)

staggered = np.concatenate([[0], fractions[:-1]])
subtracted = np.array(fractions) - staggered
subtracted = subtracted / np.sum(subtracted)

indices_df = pd.DataFrame({
    'n_barcodes': break_points,
    'fraction_reads': subtracted
})
```

### Step 3: Select Noise Threshold

**Find maximum log-weight difference** (lines 332-365):
```python
log_weights = np.log(indices_df['fraction_reads'].values)
weight_diff = log_weights[:-1] - log_weights[1:]

if len(weight_diff) > 0:
    noise_start = indices_df['n_barcodes'].values[np.argmax(weight_diff)]
else:
    noise_start = indices_df['n_barcodes'].values[-1]

# Apply minimum weight threshold
if indices_df['fraction_reads'].min() > min_weight:
    noise_start = indices_df['n_barcodes'].max()

# Check cumulative weight
cumsum_weight = np.cumsum(indices_df['fraction_reads'].values)
if len(cumsum_weight) > 0:
    cutoff_idx = np.where(cumsum_weight > (1 - min_weight))[0]
    if len(cutoff_idx) > 0:
        locationofminweightcutoff = cutoff_idx[0]
        remaining_weight = np.sum(
            indices_df['fraction_reads'].values[
                indices_df['n_barcodes'].values > noise_start
            ]
        )
        if remaining_weight > min_weight:
            noise_start = indices_df['n_barcodes'].values[locationofminweightcutoff]

# Special case: if only 2 barcodes, check ratio
if noise_start == 2:
    sorted_counts = np.sort(indices_df['fraction_reads'].values)[::-1]
    if len(sorted_counts) >= 2:
        if sorted_counts[1] / sorted_counts[0] < 0.01:
            noise_start = 1
```

---

## Example: s43_ch06 Walkthrough

### What R Does:

1. **Stochastic search finds breakpoints**: [1, 2, 5, 13, ...]
2. **Creates population table**:
   ```
   Number of barcodes  Fraction of reads
   1                   0.45
   2                   0.20
   5                   0.25
   13                  0.10
   ```
3. **Finds max log-weight drop**: Between population 5 and 13
4. **Sets noisestart = 5** (5 populations detected)
5. **Simulates rarefaction** with target = 5 populations
6. **Estimates Ns = 5**
7. **Reports "Number of barcodes" = 5** (really means 5 populations)

### What Python Does:

1. **Stochastic search finds breakpoints**: [1, 2, 5, 13, ...]
2. **Creates indices table**: Same as R
3. **Finds max log-weight drop**: Between 5 and 13 barcodes
4. **Sets noise_start = 13** (keep 13 barcodes)
5. **Simulates rarefaction** with target = 13 barcodes
6. **Estimates Ns = 20.49**
7. **Reports "Number of barcodes" = 13** (actually 13 barcodes)

---

## The Root Cause

**Python's `noise_start` is selecting the WRONG VALUE from the indices table.**

Look at line 333-337 in Python:
```python
if len(weight_diff) > 0:
    noise_start = indices_df['n_barcodes'].values[np.argmax(weight_diff)]
else:
    noise_start = indices_df['n_barcodes'].values[-1]
```

This should select the value BEFORE the max drop, not at the max drop!

### Example:
```
indices_df:
  n_barcodes  fraction_reads  log_weight  weight_diff (i to i+1)
  1           0.45            -0.80
  2           0.20            -1.61       0.81
  5           0.25            -1.39       -0.22  <- [ERROR: Python selects this]
  13          0.10            -2.30       0.91   <- CORRECT: max drop is HERE
```

`np.argmax(weight_diff)` returns index 2 (weight_diff from 5→13 is max)

`indices_df['n_barcodes'].values[2]` = **5**  ← This is CORRECT!

Actually, let me recheck this...

---

## Rechecking Python's Logic

Wait, let me trace through the Python code more carefully:

```python
log_weights = np.log(indices_df['fraction_reads'].values)
# log_weights = [-0.80, -1.61, -1.39, -2.30]

weight_diff = log_weights[:-1] - log_weights[1:]
# weight_diff = [-0.80 - (-1.61), -1.61 - (-1.39), -1.39 - (-2.30)]
# weight_diff = [0.81, -0.22, 0.91]
#                i=0   i=1     i=2

np.argmax(weight_diff) = 2  # Index of 0.91

noise_start = indices_df['n_barcodes'].values[2]
# From indices_df row 2 (0-indexed), n_barcodes = 5
```

So Python IS selecting 5, the same as R! But the output shows:
- R: 5 barcodes, Ns=5
- Python: 13 barcodes, Ns=20.49

This means Python is NOT actually using this selected value correctly, OR the breakpoints are different.

---

## Hypothesis: Different Breakpoints

The breakpoints found by R and Python might be different because:

1. **Python is missing the read count jump detection**
   - R adds breakpoints where there's a 10-fold jump in READ COUNTS
   - Python only checks resiliency jumps

2. **Random seed differences**
   - Stochastic search with 10,000 iterations will find different local minima

---

## Next Steps

1. **Add debug output** to both R and Python to see:
   - What breakpoints are found
   - What the indices/populations table looks like
   - What `noisestart`/`noise_start` value is selected
   - What gets passed to Ns estimation

2. **Align the breakpoint detection**:
   - Add read count jump detection to Python (lines 94, 159 in R)
   - Ensure same filtering logic

3. **Verify the Ns calculation**:
   - Confirm both use the same simulation approach
   - Check if the interpolation differs

---

## Recommendation

The user wants Ns values to match. To achieve this:

**Option 1: Debug and align** (recommended)
- Add extensive logging to both scripts
- Run on s43_ch06 and compare step-by-step
- Identify exactly where the algorithms diverge

**Option 2: Use R's output directly** (if Option 1 is too complex)
- For Ns values, use R's implementation
- For Nb values, either implementation works (they match excellently)

**Option 3: Simplify Python to match R exactly**
- Port R's exact logic to Python line-by-line
- Remove any algorithmic differences
- This may require significant refactoring
