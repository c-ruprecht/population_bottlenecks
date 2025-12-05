"""
Nb Calculator - Bottleneck Population Size Estimation for Barcode Sequencing Data

This script estimates the bottleneck population size (Nb) from barcode sequencing
experiments using Wright-Fisher F-statistic methods.

Tested against stampr input outputs from rtisan repo and Nb matches perfectly for lebruin 2024 dataset

Usage:
    python getFP_Nb.py <reads_table.csv> <reference_columns> <output_dir>

Example:
    python getFP_Nb.py data/readstable.csv "0" output/
"""

import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def calculate_wright_fisher_bottleneck(input_vec, output_vec):
    """
    Calculate bottleneck population size (Nb) using Wright-Fisher F-statistic.

    This estimates the effective population size from frequency changes between
    input and output populations using variance in allele frequencies.

    Parameters:
    -----------
    input_vec : array
        Reference population barcode counts
    output_vec : array
        Test sample barcode counts

    Returns:
    --------
    Nb : float
        Bottleneck population size estimate
    F : float
        F-statistic (variance ratio)
    """
    input_prop = input_vec / np.sum(input_vec)
    output_prop = output_vec / np.sum(output_vec)

    # Calculate variance components
    numerator = (output_prop - input_prop) ** 2
    denominator = input_prop * (1 - input_prop)

    # Remove infinite values
    sigma = numerator / denominator
    sigma = sigma[~np.isinf(sigma)]

    # F-statistic
    F = np.nanmean(sigma)

    # Bottleneck size estimate (accounting for sampling variance)
    Nb = 1 / (F - 1/np.sum(input_vec) - 1/np.sum(output_vec))

    return max(Nb, 1/F), F


def remove_low_frequency_noise(output_vec, min_freq=1e-7):
    """Remove extremely low-frequency barcodes that are likely noise."""
    output_freq = output_vec / np.sum(output_vec)
    output_vec_clean = output_vec.copy()
    output_vec_clean[output_freq < min_freq] = 0
    return output_vec_clean


def adjust_singleton_noise(output_vec):
    """
    Correct for singleton barcode noise using ratio-based thresholding.

    If singletons (count=1) are much more abundant than count>1 barcodes,
    subtract 1 from all counts to remove likely sequencing errors.
    """
    # If singletons dominate (ratio < 1.2), apply noise correction
    if np.sum(output_vec == 1) > 1.2 * np.sum(output_vec > 1):
        output_vec_corrected = output_vec - 1
        output_vec_corrected[output_vec_corrected < 0] = 0
        return output_vec_corrected

    return output_vec


def analyze_sample(sample_name, input_vec, output_vec):
    """
    Calculate Nb for a single sample.

    Parameters:
    -----------
    sample_name : str
        Sample identifier
    input_vec : array
        Reference population barcode counts
    output_vec : array
        Sample barcode counts

    Returns:
    --------
    results : dict
        Dictionary containing Nb and F-statistic
    """
    print(f"\nAnalyzing sample: {sample_name}")

    # Store original total reads
    total_reads = np.sum(output_vec)

    # Noise correction
    output_vec = remove_low_frequency_noise(output_vec)
    output_vec = adjust_singleton_noise(output_vec)

    # Calculate bottleneck size
    Nb, F_stat = calculate_wright_fisher_bottleneck(input_vec, output_vec)

    results = {
        'sample_name': sample_name,
        'total_reads': total_reads,
        'Nb': Nb,
        'F_statistic': F_stat,
        'n_barcodes_detected': np.sum(output_vec > 0)
    }

    print(f"  Nb = {Nb:.2f}")
    print(f"  F-statistic = {F_stat:.6f}")
    print(f"  Barcodes detected = {results['n_barcodes_detected']}")

    return results


def main():
    """Main entry point for command-line execution."""
    parser = argparse.ArgumentParser(
        description='Calculate bottleneck population size (Nb) from barcode sequencing'
    )
    parser.add_argument('reads_table', type=str,
                       help='Path to reads table CSV (barcodes × samples)')
    parser.add_argument('reference_columns', type=str,
                       help='Comma-separated column indices for reference (0-indexed)')
    parser.add_argument('output_dir', type=str,
                       help='Output directory path')

    args = parser.parse_args()

    # Load data
    reads_table = pd.read_csv(args.reads_table, index_col=0)

    # Parse reference columns
    ref_cols = [int(x) for x in args.reference_columns.split(',')]

    # Calculate reference vector
    reference_vec = reads_table.iloc[:, ref_cols].sum(axis=1).values

    # Get sample columns
    sample_cols = [i for i in range(len(reads_table.columns)) if i not in ref_cols]
    sample_names = reads_table.columns[sample_cols]

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Analyze each sample
    all_results = []

    for sample_name in sample_names:
        output_vec = reads_table[sample_name].values

        results = analyze_sample(sample_name, reference_vec, output_vec)
        all_results.append(results)

    # Create summary table
    summary_df = pd.DataFrame([
        {
            'sample': r['sample_name'],
            'total_reads': r['total_reads'],
            'n_barcodes': r['n_barcodes_detected'],
            'Nb': r['Nb'],
            'F_statistic': r['F_statistic']
        }
        for r in all_results
    ])

    # Save outputs
    summary_df.to_csv(output_dir / 'Nb_estimates.csv', index=False)

    print("\n" + "="*60)
    print("Analysis complete!")
    print(f"Results saved to: {output_dir / 'Nb_estimates.csv'}")
    print("="*60)
    print("\nSummary:")
    print(summary_df.to_string(index=False))


if __name__ == '__main__':
    main()
