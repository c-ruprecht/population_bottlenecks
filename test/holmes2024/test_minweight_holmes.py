#!/usr/bin/env python3
"""
Test different minweight values for holmes2024 dataset.

Expected R output for s43_ch06:
- Ns = 4
- Nb = 30.8689044295878
- Number of barcodes = 4
"""

import subprocess
import pandas as pd
import sys
import os

# Target values from R
TARGET_NS = 4
TARGET_NB = 30.8689044295878
TARGET_BARCODES = 4
SAMPLE_NAME = "s43_ch06"

# Test different minweight values
minweight_values = [0.99, 0.97, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01]

print("="*80)
print(f"Testing minweight values to match R output for {SAMPLE_NAME}")
print("="*80)
print(f"\nTarget values from R:")
print(f"  Ns: {TARGET_NS}")
print(f"  Nb: {TARGET_NB:.6f}")
print(f"  Number of barcodes: {TARGET_BARCODES}")
print("\n" + "="*80)

results = []

for minweight in minweight_values:
    output_dir = f"test/holmes2024/minweight_test/mw_{minweight}"
    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        "python", "scripts/getFP.py",
        "test/holmes2024/inputs/NoHopFreq_Masked_tatC_subset.csv",
        "test/holmes2024/inputs/CFU_tatC_12112023.csv",
        "0,1,2,3,4,5",  # Gavage columns (0-indexed)
        str(minweight),
        output_dir
    ]

    print(f"\nTesting minweight = {minweight}...", end=" ", flush=True)

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"ERROR")
        print(f"  STDERR: {result.stderr[:200]}")
        continue

    # Read the results
    estimates_file = os.path.join(output_dir, "TableOfEstimates.csv")

    if not os.path.exists(estimates_file):
        print("NO OUTPUT")
        continue

    df = pd.read_csv(estimates_file)

    # Find the row for our sample
    sample_row = df[df['sample'] == SAMPLE_NAME]

    if len(sample_row) == 0:
        print("SAMPLE NOT FOUND")
        continue

    ns = sample_row['Ns'].values[0]
    nb = sample_row['Nb'].values[0]
    n_barcodes = sample_row['Number of barcodes'].values[0]

    # Calculate differences
    ns_diff = abs(ns - TARGET_NS)
    nb_diff = abs(nb - TARGET_NB)
    bc_diff = abs(n_barcodes - TARGET_BARCODES)

    print(f"Ns={ns:.2f} (diff={ns_diff:.2f}), Nb={nb:.6f} (diff={nb_diff:.6f}), Barcodes={n_barcodes} (diff={bc_diff})")

    results.append({
        'minweight': minweight,
        'Ns': ns,
        'Nb': nb,
        'Barcodes': n_barcodes,
        'Ns_diff': ns_diff,
        'Nb_diff': nb_diff,
        'Barcodes_diff': bc_diff
    })

    # Check if we found a perfect match for Ns and barcodes
    if ns == TARGET_NS and n_barcodes == TARGET_BARCODES:
        print(f"  *** PERFECT MATCH for Ns and Barcodes! ***")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)

results_df = pd.DataFrame(results)

# Sort by combined difference
results_df['total_diff'] = results_df['Ns_diff'] + results_df['Barcodes_diff'] + results_df['Nb_diff']/100
results_df = results_df.sort_values('total_diff')

print(f"\n{'minweight':<12} {'Ns':<10} {'Nb':<12} {'Barcodes':<10} {'Ns_diff':<10} {'Nb_diff':<12} {'BC_diff':<10}")
print("-"*80)

for _, row in results_df.iterrows():
    marker = " ***" if row['Ns'] == TARGET_NS and row['Barcodes'] == TARGET_BARCODES else ""
    print(f"{row['minweight']:<12.2f} {row['Ns']:<10.2f} {row['Nb']:<12.6f} {row['Barcodes']:<10.0f} {row['Ns_diff']:<10.2f} {row['Nb_diff']:<12.6f} {row['Barcodes_diff']:<10.0f}{marker}")

# Find best match
best_match = results_df.iloc[0]
print(f"\n{'='*80}")
print(f"BEST MATCH:")
print(f"  minweight = {best_match['minweight']}")
print(f"  Ns = {best_match['Ns']:.2f} (R: {TARGET_NS})")
print(f"  Nb = {best_match['Nb']:.6f} (R: {TARGET_NB:.6f})")
print(f"  Barcodes = {best_match['Barcodes']:.0f} (R: {TARGET_BARCODES})")
print(f"{'='*80}")

# Save results
results_df.to_csv('test/holmes2024/minweight_test/minweight_comparison.csv', index=False)
print(f"\nResults saved to: test/holmes2024/minweight_test/minweight_comparison.csv")
