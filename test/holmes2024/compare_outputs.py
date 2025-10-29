#!/usr/bin/env python3
"""Compare R and Python getFP outputs for holmes2024 dataset."""

import pandas as pd
import numpy as np

# Read R and Python outputs
r_file = 'test/holmes2024/outputs/NsNb_tatC_12112023.csv'
py_file = 'test/holmes2024/python_outputs/TableOfEstimates.csv'

r_df = pd.read_csv(r_file, index_col=0)
py_df = pd.read_csv(py_file)

# Clean sample names
r_df.index = r_df.index.str.strip('"')
py_df['sample'] = py_df['sample'].str.strip()

print("="*80)
print("HOLMES2024 DATASET: R vs Python getFP Comparison")
print("="*80)

print(f"\nR output: {len(r_df)} samples")
print(f"Python output: {len(py_df)} samples")

# Align samples
common_samples = sorted(set(r_df.index) & set(py_df['sample']))
print(f"Common samples: {len(common_samples)}")

print("\n" + "="*80)
print("SAMPLE-BY-SAMPLE COMPARISON")
print("="*80)

print(f"\n{'Sample':<15} {'R_Nb':<12} {'Py_Nb':<12} {'Nb_Diff':<12} {'R_Ns':<10} {'Py_Ns':<10} {'Ns_Diff':<10}")
print("-"*85)

nb_diffs = []
ns_diffs = []
barcode_diffs = []

for sample in common_samples:
    r_row = r_df.loc[sample]
    py_row = py_df[py_df['sample'] == sample].iloc[0]

    r_nb = r_row['Nb']
    py_nb = py_row['Nb']
    r_ns = r_row['Ns']
    py_ns = py_row['Ns']

    nb_diff = abs(r_nb - py_nb)
    ns_diff = abs(r_ns - py_ns)

    nb_diffs.append(nb_diff)
    ns_diffs.append(ns_diff)

    print(f"{sample:<15} {r_nb:<12.4f} {py_nb:<12.4f} {nb_diff:<12.4f} {r_ns:<10.2f} {py_ns:<10.2f} {ns_diff:<10.2f}")

print("\n" + "="*80)
print("STATISTICAL SUMMARY")
print("="*80)

print(f"\nNb (Bottleneck Size) Comparison:")
print(f"  Mean absolute difference: {np.mean(nb_diffs):.6f}")
print(f"  Median absolute difference: {np.median(nb_diffs):.6f}")
print(f"  Max absolute difference: {np.max(nb_diffs):.6f}")
print(f"  Min absolute difference: {np.min(nb_diffs):.6f}")

# Percent differences
r_nb_values = [r_df.loc[s]['Nb'] for s in common_samples]
nb_percent_diffs = [100 * abs(nb_diffs[i] / r_nb_values[i]) if r_nb_values[i] != 0 else 0
                    for i in range(len(nb_diffs))]
print(f"  Mean percent difference: {np.mean(nb_percent_diffs):.2f}%")
print(f"  Median percent difference: {np.median(nb_percent_diffs):.2f}%")

# Matches
close_matches_nb = sum(1 for d in nb_diffs if d < 1)
very_close_matches_nb = sum(1 for d in nb_diffs if d < 0.1)
print(f"  Samples with Nb diff < 1: {close_matches_nb}/{len(common_samples)} ({100*close_matches_nb/len(common_samples):.1f}%)")
print(f"  Samples with Nb diff < 0.1: {very_close_matches_nb}/{len(common_samples)} ({100*very_close_matches_nb/len(common_samples):.1f}%)")

print(f"\nNs (Founding Population) Comparison:")
print(f"  Mean absolute difference: {np.mean(ns_diffs):.2f}")
print(f"  Median absolute difference: {np.median(ns_diffs):.2f}")
print(f"  Max absolute difference: {np.max(ns_diffs):.2f}")
print(f"  Min absolute difference: {np.min(ns_diffs):.2f}")

# Percent differences for Ns
r_ns_values = [r_df.loc[s]['Ns'] for s in common_samples]
ns_percent_diffs = [100 * abs(ns_diffs[i] / r_ns_values[i]) if r_ns_values[i] != 0 else 0
                    for i in range(len(ns_diffs))]
print(f"  Mean percent difference: {np.mean(ns_percent_diffs):.2f}%")
print(f"  Median percent difference: {np.median(ns_percent_diffs):.2f}%")

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)

if np.mean(nb_percent_diffs) < 5:
    print("\n✓✓✓ Nb values match EXCELLENTLY (mean diff < 5%)")
elif np.mean(nb_percent_diffs) < 10:
    print("\n✓✓ Nb values match VERY WELL (mean diff < 10%)")
elif np.mean(nb_percent_diffs) < 20:
    print("\n✓ Nb values match REASONABLY (mean diff < 20%)")
else:
    print("\n✗ Nb values show significant differences")

if np.mean(ns_percent_diffs) < 10:
    print("✓✓✓ Ns values match EXCELLENTLY (mean diff < 10%)")
elif np.mean(ns_percent_diffs) < 50:
    print("✓ Ns values match REASONABLY (mean diff < 50%)")
else:
    print("✗ Ns values show significant differences (expected due to algorithm differences)")

print("\n" + "="*80)
print("\nNote: Nb (bottleneck size) is the key metric for population bottleneck analysis.")
print("Ns differences are expected due to different population detection algorithms.")
print("="*80)

# Save comparison
comparison_df = pd.DataFrame({
    'sample': common_samples,
    'R_Nb': [r_df.loc[s]['Nb'] for s in common_samples],
    'Python_Nb': [py_df[py_df['sample'] == s].iloc[0]['Nb'] for s in common_samples],
    'Nb_diff': nb_diffs,
    'Nb_percent_diff': nb_percent_diffs,
    'R_Ns': [r_df.loc[s]['Ns'] for s in common_samples],
    'Python_Ns': [py_df[py_df['sample'] == s].iloc[0]['Ns'] for s in common_samples],
    'Ns_diff': ns_diffs,
    'Ns_percent_diff': ns_percent_diffs
})

comparison_df.to_csv('test/holmes2024/python_outputs/comparison_R_vs_Python.csv', index=False)
print(f"\nComparison saved to: test/holmes2024/python_outputs/comparison_R_vs_Python.csv")
