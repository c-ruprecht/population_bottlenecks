# MinWeight Analysis: Why Ns Values Differ Between R and Python

## Question
What minweight parameter was used in the R script to produce the paper outputs?

## Experiment
Created a subset of ml01_NoHopping.csv with:
- All 26 gavage/input samples (columns 0-26)
- One test sample: M1024_S11_L001_R1_001

Tested minweight values: 0.97, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01

## Results

### Target (R Output for M1024_S11_L001_R1_001):
- **Ns**: 3
- **Nb**: 6.739492253
- **Number of barcodes**: 3

### Python Output (ALL minweight values):
- **Ns**: 40.82
- **Nb**: 6.739492
- **Number of barcodes**: 44

### Key Observations:

1. **Nb matches PERFECTLY** (6.739492) across ALL minweight values
   - This confirms the bottleneck size calculation is correct
   - The core genetic distance algorithm is working properly

2. **Ns and barcode count NEVER change** regardless of minweight
   - Python always gives: Ns=40.82, Barcodes=44
   - R gives: Ns=3, Barcodes=3
   - Difference: Python keeps 14.7x more barcodes than R

3. **minweight parameter is NOT the source of the difference**
   - Varying minweight from 0.01 to 0.97 had ZERO effect
   - The difference must be in a different part of the algorithm

## Investigation of R's FrequenciesWithoutNoise Output

Checked the R paper output file `FrequenciesWithoutNoise2023.csv`:
- **Non-zero barcodes for M1024**: 201 barcodes

But the NsNb table says:
- **Number of barcodes**: 3

### This reveals a KEY insight:

The "Number of barcodes" in the NsNb table is NOT the count of barcodes in FrequenciesWithoutNoise!

**Two possibilities:**

1. **FrequenciesWithoutNoise is an intermediate output**
   - The getFP.R function applies additional filtering AFTER creating this file
   - The final barcode count (3) comes from more aggressive population detection
   - The "Number of barcodes" represents detected populations, not total barcodes

2. **Different definition of "barcode"**
   - FrequenciesWithoutNoise shows all non-noise barcodes (201)
   - NsNb "Number of barcodes" might be the number of FOUNDER barcodes
   - Or the number detected at a specific population threshold

## Conclusion

### Finding:
**The minweight parameter is NOT causing the Ns difference.**

The difference between R and Python is in:
1. **Population detection logic** - How many distinct populations are identified
2. **Founder population estimation** - How Ns is calculated from the detected populations
3. **The "Number of barcodes" metric** - This appears to be a derived metric, not a simple count

### Why Nb Matches Perfectly:
- Nb (bottleneck size) is calculated from the initial genetic distance
- This calculation is identical between R and Python
- The noise filtering doesn't affect Nb calculation

### Why Ns Differs:
- Ns (founding population size) depends on how populations are detected and counted
- R's getFP appears to use more complex population detection
- The Python implementation may be using a simpler counting method
- This is likely in the `determine_noise_threshold()` or resiliency curve analysis

## Recommendation

To understand the true difference, need to:
1. **Examine the R getFP.R script's population detection logic more carefully**
   - Look at how it determines "populations" vs "barcodes"
   - Check if there's post-processing after FrequenciesWithoutNoise is written

2. **Add debug output to both R and Python**
   - Print intermediate values during population detection
   - Track how many populations vs barcodes are identified at each step

3. **Compare the resiliency curve calculation**
   - This may be where populations are being merged or filtered

## Status
- ✓ minweight is NOT the parameter causing differences
- ✗ The actual difference is in population detection/counting logic
- → Need deeper investigation of the algorithm logic, not parameter tuning
