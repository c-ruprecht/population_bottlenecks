# R vs Python Ns Comparison for holmes2024

## Environment Setup

Created conda environment `r_python_env` with both R 4.2 and Python 3.9 installed.

Location: `/Users/ruprec01/opt/anaconda3/envs/r_python_env`

To activate:
```bash
conda activate r_python_env
```

##Results Summary

### Sample: s43_ch06

#### R Paper Output (from NsNb_tatC_12112023.csv):
- **Ns**: 4
- **Nb**: 30.8689044295878
- **Number of barcodes**: 4
- **TotalReads**: 772,382

#### Python getFP_fixNs.py Output:
- **Ns**: 20.49
- **Nb**: 31.173166
- **Number of barcodes**: 13
- **TotalReads**: 772,382

### Comparison:

| Metric | R | Python | Difference | % Difference |
|--------|---|--------|------------|--------------|
| **Nb** | 30.869 | 31.173 | 0.304 | 0.99% |
| **Ns** | 4 | 20.49 | 16.49 | 412% |
| **Barcodes** | 4 | 13 | 9 | 225% |

## Analysis

### 1. Nb Matches Excellently ✓
- Difference: 0.304 (0.99%)
- This confirms the genetic distance calculation is correct
- Bottleneck size estimation is working properly

### 2. Ns and Barcode Count Differ Significantly
- Python keeps 13 barcodes vs R's 4
- Python estimates Ns=20.49 vs R's Ns=4
- This suggests much more aggressive noise filtering in R

### 3. What "Number of Barcodes" Means

This is the KEY question. Two possibilities:

**Hypothesis 1**: Different noise filtering thresholds
- R filters out more low-abundance barcodes
- Python keeps barcodes that R considers "noise"
- Result: Different barcode counts

**Hypothesis 2**: "Number of barcodes" = "Number of populations"
- R's "Number of barcodes" might be counting detected **populations**, not total barcodes
- The R FrequenciesWithoutNoise file might contain ALL noise-filtered barcodes
- The NsNb table's "Number of barcodes" is a derived metric from population detection
- This would explain why the lebruin2024 dataset showed 201 barcodes in FrequenciesWithoutNoise but only 3 in NsNb

## Next Steps to Match R Output

### Option 1: Examine R's FrequenciesWithoutNoise Output
Compare the R paper's FrequenciesWithoutNoise_tatC_12112023.csv to understand:
- How many barcodes does R keep for s43_ch06 after noise filtering?
- Does this match the "4 barcodes" in NsNb table?

### Option 2: Add Debug Output to Both Scripts
Modify both scripts to print:
- Number of barcodes after each filtering step
- Population detection breakpoints
- Which barcodes are kept vs filtered

### Option 3: Study the R Code More Carefully
Focus on these sections in getFP_12112023.R:
1. **Noise filtering logic** (where barcodes are removed)
2. **Population detection** (how distinct populations are identified)
3. **Ns calculation** (how founding population is estimated)
4. **Relationship between barcodes and Ns**

## Current Understanding

The difference is NOT in:
- ✓ minweight parameter (tested 0.01-0.99, no effect)
- ✓ Nb calculation (matches within 1%)
- ✓ Core genetic distance algorithm (GD matches 93%)

The difference IS in:
- ✗ Noise filtering aggressiveness
- ✗ Population structure detection
- ✗ Interpretation of "Number of barcodes"

## Recommendation for Moving Forward

To match R's Ns values, we need to understand:

1. **What does R's "Number of barcodes" actually represent?**
   - Total barcodes after noise filtering?
   - Number of detected populations?
   - Number of founder barcodes?

2. **How does R go from barcode counts to Ns?**
   - Is it a direct count?
   - Is it estimated from resiliency curves?
   - Is there population merging/detection?

3. **Check R's intermediate outputs**
   - FrequenciesWithoutNoise for s43_ch06
   - Resiliency curve calculations
   - Population breakpoints

## Files for Comparison

**R Outputs (Paper)**:
- `test/holmes2024/outputs/NsNb_tatC_12112023.csv`
- `test/holmes2024/outputs/FrequenciesWithoutNoise_tatC_12112023.csv`

**Python Outputs**:
- `test/holmes2024/python_fixNs_output/TableOfEstimates.csv`
- `test/holmes2024/python_outputs/TableOfEstimates.csv`

**Scripts**:
- `test/holmes2024/getFP_12112023.R`
- `test/holmes2024/getFP_fixNs.py` (same as scripts/getFP.py)

## Action Items

1. ✓ Create R+Python environment
2. ✗ Run R script on subset (failed - needs debugging)
3. ⚠️ Compare R's FrequenciesWithoutNoise to understand barcode filtering
4. ⚠️ Add detailed logging to both R and Python scripts
5. ⚠️ Identify exact line where R and Python diverge in their calculations
