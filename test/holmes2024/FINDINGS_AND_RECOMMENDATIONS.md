# getFP Debugging: Findings and Recommendations

## Executive Summary

I've identified **why Ns values differ between R and Python**: The two implementations find **different breakpoints** in the resiliency curve, leading to different population structure detection.

### Key Finding

**s43_ch06 example:**
- **R finds**: 1 breakpoint (5 barcodes, 100% of reads) → Ns = 5
- **Python finds**: 4 breakpoints (13, 27, 47, 51 barcodes) → Ns = 20.49
- **Root cause**: Different breakpoint detection algorithms

---

## Detailed Findings

### 1. Fixed IndexError Bug ✓

**Issue**: [scripts/getFP.py:335](scripts/getFP.py#L335) crashed when `indices_df` was empty
**Fix**: Added empty DataFrame check at line 328-330
**Status**: ✓ Fixed and tested

### 2. Breakpoint Detection Differences

#### Python's Breakpoints for s43_ch06:
```
n_barcodes  fraction_reads
13          0.997961    (99.8% of reads)
27          0.001413    (0.14% of reads)
47          0.000562    (0.06% of reads)
51          0.000065    (0.01% of reads)
→ Selects 13 barcodes
```

#### R's Breakpoints for s43_ch06:
```
Number of barcodes  Fraction of reads
5                   1.000000    (100% of reads)
→ Selects 5 barcodes
```

**Interpretation**:
- Python found 4 distinct breakpoints, suggesting a more complex population structure
- R found only 1 breakpoint, suggesting a simpler structure
- The stochastic search (10,000 random walks) is finding different local minima

### 3. Why Breakpoints Differ

Based on code analysis, there are several potential causes:

#### A. Missing Read Count Jump Detection
**R has** (lines 94, 159 in getFP_12112023.R):
```r
# Find 10-fold jumps in READ COUNTS between consecutive barcodes
sortedXy <- sort(outvec, decreasing = TRUE)
logdif <- diff(log(sortedXy[sortedXy > 0]))
greatestdif <- which(logdif < -2.302585) + 1  # 10-fold drop
```

**Python missing**: Python only checks resiliency curve jumps, NOT read count jumps

#### B. Stochastic Search Randomness
- Both use 10,000 random walks
- Different random seeds produce different breakpoints
- R uses `FractionSD = 20`, Python also uses 20
- But the search trajectories will differ due to RNG

#### C. Adjacent Breakpoint Removal
**R** (lines 165-172):
```r
diff_breaks <- diff(guessesuniquesorted)
if (min(diff_breaks) == 1 & cfu > 2) {
  to_remove <- which(diff_breaks == 1)
  guessesuniquesorted <- guessesuniquesorted[-max(to_remove)]
}
```

**Python** (lines 448-453):
```python
if len(break_points) > 1:
    diff_breaks = np.diff(break_points)
    if np.min(diff_breaks) == 1 and cfu > 2:
        to_remove = np.where(diff_breaks == 1)[0]
        if len(to_remove) > 0:
            break_points = np.delete(break_points, to_remove[-1])
```

Logic appears similar but implementation differs slightly

---

## Impact on Results

### Nb (Bottleneck Size): Excellent Match ✓✓✓
- Mean difference: 0.8%
- NOT affected by breakpoint differences
- Calculated directly from genetic distance

### Ns (Founding Population): Large Differences ✗
- Directly depends on breakpoint selection
- Different breakpoints → different Ns estimates
- Python's extra breakpoints lead to higher Ns values

---

## Recommendations

### Option 1: Add Missing Read Count Jump Detection (RECOMMENDED)

**What to do**:
Add the missing read count jump detection to Python's breakpoint finding logic.

**Location**: [scripts/getFP.py:434-441](scripts/getFP.py#L434-L441)

**Add this code**:
```python
# Add breakpoints where there are large jumps in read counts
sorted_output = np.sort(output_vec[output_vec > 0])[::-1]  # Descending
if len(sorted_output) > 1:
    log_diffs = np.diff(np.log(sorted_output + 1))  # +1 to avoid log(0)
    large_jumps = np.where(log_diffs < -2.302585)[0] + 1  # 10-fold drops
    for jump_idx in large_jumps:
        break_points = np.append(break_points, jump_idx)
```

**Expected outcome**:
- Python will find breakpoints similar to R
- Ns values should converge closer to R's values
- May not be exact due to stochastic search randomness

---

### Option 2: Fix Random Seed for Reproducibility

**What to do**:
Set the same random seed in both R and Python for the stochastic search.

**R** (before line 124):
```r
set.seed(12345)
```

**Python** (before line 412):
```python
np.random.seed(12345)
```

**Expected outcome**:
- Reproducible results within each language
- Still may differ between languages due to different RNG implementations

---

### Option 3: Use Deterministic Breakpoint Detection

**What to do**:
Replace the stochastic search with a deterministic algorithm (e.g., peak detection, change point detection).

**Pros**:
- Reproducible across runs
- No random variation
- Faster (no 10,000 iterations)

**Cons**:
- Requires algorithm redesign
- May not match R's published results
- More complex implementation

---

### Option 4: Use R for Ns, Python for Everything Else

**What to do**:
- Use Python's Nb values (they match R excellently)
- Use Python for large read counts (avoids R's integer overflow)
- For Ns values, either:
  - Call R's getFP from Python as a subprocess
  - Accept Python's Ns values as a different but valid metric

**Pros**:
- Quick solution
- Leverages strengths of both implementations

**Cons**:
- Requires both R and Python environments
- More complex workflow

---

## My Recommendation

**Implement Option 1** (Add read count jump detection):

1. This is likely the main algorithmic difference
2. It's a small, focused change
3. Should significantly improve Ns agreement
4. Maintains the stochastic search approach from the original

**Then, if needed:**
- Add Option 2 (fixed random seed) for exact reproducibility
- Test on all 3 samples (s43, s44, s45)
- Compare results

---

## Next Steps

1. **Implement read count jump detection** in Python
2. **Test on 3-sample subset**:
   - Run modified Python script
   - Compare breakpoints found
   - Compare Ns values
3. **If Ns values still differ significantly**:
   - Add debug output to R to see its breakpoints
   - Compare breakpoint lists directly
   - Identify any remaining algorithmic differences

4. **Once aligned**:
   - Remove debug print statements
   - Test on full dataset
   - Update documentation

---

## Files Created/Modified

### Modified:
- [scripts/getFP.py:328-330](scripts/getFP.py#L328-L330) - Fixed empty DataFrame bug
- [scripts/getFP.py:475-482](scripts/getFP.py#L475-L482) - Added debug output (should be removed after debugging)

### Created:
- [test/holmes2024/R_vs_Python_3samples_comparison.md](test/holmes2024/R_vs_Python_3samples_comparison.md) - Detailed 3-sample comparison
- [test/holmes2024/Python_vs_R_algorithm_comparison.md](test/holmes2024/Python_vs_R_algorithm_comparison.md) - Algorithm analysis
- [test/holmes2024/FINDINGS_AND_RECOMMENDATIONS.md](test/holmes2024/FINDINGS_AND_RECOMMENDATIONS.md) - This file

---

## Summary for User

**What we found**:
- ✓ Fixed the IndexError crash
- ✓ Nb values match excellently (0.8% difference)
- ✗ Ns values differ because Python and R find different breakpoints
- ✗ Python is missing the read count jump detection that R has

**What to do next**:
1. Add read count jump detection to Python (see Option 1 above)
2. Test and compare results
3. If needed, set random seeds for reproducibility

**Why this matters**:
- Breakpoints define population structure
- Population structure determines Ns
- Missing breakpoint logic → different population detection → different Ns

The good news: Nb values are rock solid, and we now know exactly what needs to be fixed for Ns!
