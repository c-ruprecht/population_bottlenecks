"""
getFP_clean - Simplified Founder Population Analysis for Pre-Cleaned Data

This script estimates founding population size (Ns) and bottleneck population (Nb)
from PRE-CLEANED barcode frequencies (e.g., after spike-in normalization and
background subtraction).

Use this when you've already performed:
- Spike-in normalization
- Background/control subtraction
- Quality filtering

Usage:
    python getFP_clean.py <freq_table.csv> <cfu_table.csv|NULL> <reference_columns> <output_dir>

Example:
    python getFP_clean.py data/cleaned_frequencies.csv data/cfu.csv "0" output/
"""

import numpy as np
import pandas as pd
import argparse
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def geometric_mean(data):
    """Calculate geometric mean of positive values."""
    data_positive = data[data > 0]
    if len(data_positive) == 0:
        return np.nan
    return np.exp(np.mean(np.log(data_positive)))


def calculate_wright_fisher_bottleneck(input_vec, output_vec):
    """
    Calculate bottleneck population size (Nb) using Wright-Fisher F-statistic.

    Based on allele frequency variance between input and output populations.

    Parameters:
    -----------
    input_vec : array
        Reference population barcode counts/frequencies
    output_vec : array
        Sample barcode counts/frequencies

    Returns:
    --------
    Nb : float
        Bottleneck population size estimate
    F : float
        F-statistic (standardized variance)
    """
    # Normalize to frequencies if needed
    input_freq = input_vec / np.sum(input_vec)
    output_freq = output_vec / np.sum(output_vec)

    # Wright-Fisher F-statistic
    # Variance in frequency change relative to expected binomial variance
    numerator = (output_freq - input_freq) ** 2
    denominator = input_freq * (1 - input_freq)

    # Remove infinite/invalid values
    sigma = numerator / denominator
    sigma = sigma[~np.isinf(sigma) & ~np.isnan(sigma)]

    if len(sigma) == 0:
        return np.nan, np.nan

    # F-statistic (mean standardized variance)
    F = np.mean(sigma)

    # Bottleneck size (correcting for sampling variance)
    total_input = np.sum(input_vec)
    total_output = np.sum(output_vec)

    Nb = 1 / (F - 1/total_input - 1/total_output)

    # Return max of Nb and 1/F (handles edge cases)
    return max(Nb, 1/F), F


def estimate_founding_population_rarefaction(input_vec, output_vec, n_detected_barcodes):
    """
    Estimate founding population size (Ns) using rarefaction curve.

    Simulates expected barcode detection at different founding population sizes
    and interpolates to match observed barcode richness.

    Parameters:
    -----------
    input_vec : array
        Reference population counts
    output_vec : array
        Observed sample counts
    n_detected_barcodes : int
        Number of detected barcodes in sample (non-zero)

    Returns:
    --------
    Ns : float
        Estimated founding population size
    rarefaction_curve : tuple
        (founding_sizes, expected_detections) for plotting
    """
    # Create expected distribution based on output, weighted by reference
    total_counts = np.sum(output_vec)
    reference_freq = input_vec / np.sum(input_vec)

    # First resample: multinomial based on reference frequencies
    expected_counts = np.random.multinomial(int(total_counts), reference_freq)

    # Rarefaction: sample at increasing depths
    n_unique_input = np.sum(input_vec > 0)

    # Sample from 1 to 20x the number of input barcodes
    depths_low = np.linspace(1, n_unique_input, 50, dtype=int)
    depths_high = np.linspace(n_unique_input, n_unique_input * 20, 50, dtype=int)
    founding_sizes = np.concatenate([depths_low, depths_high])

    expected_detections = []

    for n_founders in founding_sizes:
        # Simulate: how many barcodes detected with n_founders?
        # Run multiple simulations and average
        detections = []
        for _ in range(5):
            # Draw n_founders from expected distribution
            sampled = np.random.multinomial(
                n_founders,
                expected_counts / np.sum(expected_counts)
            )
            # Count how many barcodes are present (>0)
            n_present = np.sum(sampled > 0)
            detections.append(n_present)

        expected_detections.append(np.mean(detections))

    # Interpolate to find Ns that gives observed richness
    interp_func = interp1d(
        founding_sizes,
        expected_detections,
        kind='linear',
        fill_value='extrapolate'
    )

    # Fine-grained search
    search_range = np.linspace(
        np.min(founding_sizes),
        np.max(founding_sizes),
        1000
    )
    predicted_richness = interp_func(search_range)

    # Find founding size closest to observed richness
    closest_idx = np.argmin(np.abs(predicted_richness - n_detected_barcodes))
    Ns = search_range[closest_idx]

    return Ns, (founding_sizes, expected_detections)


def calculate_diversity_metrics(output_vec):
    """
    Calculate diversity metrics from barcode distribution.

    Parameters:
    -----------
    output_vec : array
        Barcode counts/frequencies

    Returns:
    --------
    metrics : dict
        Dictionary containing:
        - richness: Number of detected barcodes
        - shannon: Shannon diversity index
        - simpson: Simpson diversity index (1-D)
        - evenness: Pielou's evenness (J)
        - effective_n: Effective number of barcodes (exp(H))
        - geometric_mean_freq: Geometric mean frequency
    """
    # Filter to non-zero
    counts = output_vec[output_vec > 0]

    if len(counts) == 0:
        return {
            'richness': 0,
            'shannon': 0,
            'simpson': 0,
            'evenness': 0,
            'effective_n': 0,
            'geometric_mean_freq': np.nan
        }

    # Convert to frequencies
    freqs = counts / np.sum(counts)

    # Richness
    richness = len(counts)

    # Shannon diversity: H = -sum(p * log(p))
    shannon = -np.sum(freqs * np.log(freqs))

    # Simpson diversity: D = 1 - sum(p^2)
    simpson = 1 - np.sum(freqs ** 2)

    # Evenness: J = H / log(S)
    evenness = shannon / np.log(richness) if richness > 1 else 1

    # Effective number of species: exp(H)
    effective_n = np.exp(shannon)

    # Geometric mean frequency (for average frequency calculation)
    geom_mean_freq = geometric_mean(freqs)

    return {
        'richness': richness,
        'shannon': shannon,
        'simpson': simpson,
        'evenness': evenness,
        'effective_n': effective_n,
        'geometric_mean_freq': geom_mean_freq
    }


def analyze_sample_clean(sample_name, input_vec, output_vec, cfu=None,
                         output_dir='.', make_plots=True):
    """
    Analyze pre-cleaned barcode frequencies to estimate Nb and Ns.

    This simplified version assumes:
    - Data is already normalized (e.g., spike-in corrected)
    - Background/noise has been subtracted
    - Only real barcodes remain (no noise filtering needed)

    Parameters:
    -----------
    sample_name : str
        Sample identifier
    input_vec : array
        Reference population barcode counts
    output_vec : array
        Sample barcode counts (pre-cleaned)
    cfu : float, optional
        Colony forming units (if available)
    output_dir : str
        Output directory path
    make_plots : bool
        Whether to generate diagnostic plots

    Returns:
    --------
    results : dict
        Dictionary containing:
        - Nb: Bottleneck population size
        - Ns: Founding population size
        - F_stat: Wright-Fisher F-statistic
        - Diversity metrics (richness, Shannon, Simpson, etc.)
        - CFU/Ns ratio (if CFU provided)
    """
    print(f"\nAnalyzing sample: {sample_name}")

    # Filter to barcodes present in both input and output
    present_in_both = (input_vec > 0) & (output_vec > 0)
    input_filtered = input_vec[present_in_both]
    output_filtered = output_vec[present_in_both]

    if len(output_filtered) == 0:
        print(f"  WARNING: No shared barcodes between input and {sample_name}")
        return None

    # Calculate bottleneck size (Nb)
    Nb, F_stat = calculate_wright_fisher_bottleneck(input_filtered, output_filtered)

    # Calculate diversity metrics
    diversity = calculate_diversity_metrics(output_vec)
    n_detected = diversity['richness']

    # Calculate founding population size (Ns)
    Ns, rarefaction_curve = estimate_founding_population_rarefaction(
        input_vec, output_vec, n_detected
    )

    # Average frequency (inverse of geometric mean frequency)
    if diversity['geometric_mean_freq'] > 0:
        avg_frequency = 1 / diversity['geometric_mean_freq']
    else:
        avg_frequency = np.nan

    # Handle CFU
    cfu_value = cfu if cfu is not None else 0

    # Calculate transmission efficiency
    transmission_efficiency = (cfu_value / Ns) if (Ns > 0 and cfu_value > 0) else 0

    # Prepare results
    results = {
        'sample_name': sample_name,
        'Nb': Nb,
        'Ns': Ns,
        'F_stat': F_stat,
        'richness': diversity['richness'],
        'shannon_diversity': diversity['shannon'],
        'simpson_diversity': diversity['simpson'],
        'evenness': diversity['evenness'],
        'effective_n': diversity['effective_n'],
        'avg_frequency': avg_frequency,
        'CFU': cfu_value,
        'log10_Ns': np.log10(Ns) if Ns > 0 else 0,
        'log10_CFU': np.log10(cfu_value) if cfu_value > 0 else 0,
        'CFU_over_Ns': transmission_efficiency,
        'total_counts': np.sum(output_vec),
        'rarefaction_curve': rarefaction_curve
    }

    # Print results
    print(f"Results for {sample_name}:")
    print(f"  Nb (Bottleneck) = {Nb:.2f}")
    print(f"  Ns (Founders) = {Ns:.2f}")
    print(f"  F-statistic = {F_stat:.4f}")
    print(f"  Richness = {diversity['richness']}")
    print(f"  Shannon H' = {diversity['shannon']:.3f}")
    print(f"  Evenness = {diversity['evenness']:.3f}")
    if cfu_value > 0:
        print(f"  CFU = {cfu_value:.2e}")
        print(f"  CFU/Ns = {transmission_efficiency:.3f}")

    # Create plots if requested
    if make_plots:
        create_diagnostic_plots_clean(
            sample_name, input_vec, output_vec,
            results, output_dir
        )

    return results


def create_diagnostic_plots_clean(sample_name, input_vec, output_vec,
                                  results, output_dir):
    """
    Create diagnostic plots for cleaned frequency analysis.

    Includes:
    1. Frequency comparison (input vs output)
    2. Rarefaction curve with Ns estimate
    3. Summary statistics
    """
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # Plot 1: Input vs Output frequencies (log-log)
    ax1 = fig.add_subplot(gs[0, :])

    input_freq = input_vec / np.sum(input_vec)
    output_freq = output_vec / np.sum(output_vec)

    # Only plot barcodes present in both
    present_both = (input_vec > 0) & (output_vec > 0)

    ax1.scatter(input_freq[present_both], output_freq[present_both],
                alpha=0.5, s=30, edgecolors='black', linewidths=0.5)
    ax1.plot([1e-6, 1], [1e-6, 1], 'k--', alpha=0.3, label='No change')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Input Frequency', fontsize=12)
    ax1.set_ylabel('Output Frequency', fontsize=12)
    ax1.set_title(f'{sample_name} - Frequency Shift (Nb = {results["Nb"]:.1f})',
                  fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Rank abundance curve
    ax2 = fig.add_subplot(gs[1, 0])

    sorted_output = np.sort(output_vec[output_vec > 0])[::-1]
    ranks = np.arange(1, len(sorted_output) + 1)

    ax2.plot(ranks, sorted_output, 'o-', markersize=4, linewidth=1)
    ax2.set_yscale('log')
    ax2.set_xlabel('Rank', fontsize=11)
    ax2.set_ylabel('Abundance', fontsize=11)
    ax2.set_title('Rank Abundance Curve', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # Add diversity annotations
    ax2.text(0.95, 0.95,
             f"Richness: {results['richness']}\n"
             f"Shannon H': {results['shannon_diversity']:.2f}\n"
             f"Evenness: {results['evenness']:.2f}",
             transform=ax2.transAxes,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             fontsize=9)

    # Plot 3: Rarefaction curve
    ax3 = fig.add_subplot(gs[1, 1])

    founding_sizes, expected_detections = results['rarefaction_curve']

    ax3.plot(founding_sizes, expected_detections, 'b-', linewidth=2,
             label='Rarefaction curve')
    ax3.axhline(results['richness'], color='red', linestyle='--',
                label=f'Observed richness = {results["richness"]}')
    ax3.axvline(results['Ns'], color='green', linestyle='--',
                label=f'Ns = {results["Ns"]:.1f}')

    ax3.set_xlabel('Founding Population Size', fontsize=11)
    ax3.set_ylabel('Expected Barcodes Detected', fontsize=11)
    ax3.set_title('Rarefaction Analysis', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)

    # Plot 4: Statistics summary
    ax4 = fig.add_subplot(gs[2, :])
    ax4.axis('off')

    stats_text = f"""
    ╔══════════════════════════════════════════════════════════════════╗
    ║                    FOUNDER POPULATION ANALYSIS                    ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║  Sample: {sample_name:<54} ║
    ║                                                                  ║
    ║  BOTTLENECK ESTIMATES:                                           ║
    ║    Nb (Bottleneck size)        : {results['Nb']:>10.2f}                      ║
    ║    Ns (Founding population)    : {results['Ns']:>10.2f}                      ║
    ║    F-statistic                 : {results['F_stat']:>10.4f}                      ║
    ║                                                                  ║
    ║  DIVERSITY METRICS:                                              ║
    ║    Richness (# barcodes)       : {results['richness']:>10}                        ║
    ║    Shannon diversity (H')      : {results['shannon_diversity']:>10.3f}                      ║
    ║    Simpson diversity           : {results['simpson_diversity']:>10.3f}                      ║
    ║    Evenness (J)                : {results['evenness']:>10.3f}                      ║
    ║    Effective N (exp(H'))       : {results['effective_n']:>10.2f}                      ║
    ║                                                                  ║
    ║  TRANSMISSION:                                                   ║
    ║    CFU (if provided)           : {results['CFU']:>10.2e}                      ║
    ║    CFU/Ns ratio                : {results['CFU_over_Ns']:>10.3f}                      ║
    ║                                                                  ║
    ║  Total counts: {results['total_counts']:<47.0f} ║
    ╚══════════════════════════════════════════════════════════════════╝
    """

    ax4.text(0.05, 0.5, stats_text, fontsize=10,
            verticalalignment='center',
            family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    # Save figure
    output_path = Path(output_dir) / f'{sample_name}_FP_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved plot to {output_path}")


def main():
    """Main entry point for command-line execution."""
    parser = argparse.ArgumentParser(
        description='FP Analysis for Pre-Cleaned Barcode Frequencies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # With CFU data:
  python getFP_clean.py cleaned_freqs.csv cfu_table.csv "0" output/

  # Without CFU data:
  python getFP_clean.py cleaned_freqs.csv NULL "0,1" output/

  # Multiple reference columns:
  python getFP_clean.py cleaned_freqs.csv NULL "0,1,2" output/

Input format:
  - Rows = barcodes
  - Columns = samples
  - Reference columns specified by index (0-based)
  - Can be counts or already normalized frequencies
        """
    )

    parser.add_argument('freq_table', type=str,
                       help='Path to cleaned frequency/count table CSV')
    parser.add_argument('cfu_table', type=str,
                       help='Path to CFU table CSV or "NULL"')
    parser.add_argument('reference_columns', type=str,
                       help='Comma-separated reference column indices (0-indexed)')
    parser.add_argument('output_dir', type=str,
                       help='Output directory path')
    parser.add_argument('--no-plots', action='store_true',
                       help='Disable plot generation')

    args = parser.parse_args()

    # Load frequency table
    freq_table = pd.read_csv(args.freq_table, index_col=0)
    print(f"\nLoaded frequency table: {freq_table.shape[0]} barcodes × {freq_table.shape[1]} samples")

    # Load CFU table if provided
    cfu_dict = {}
    if args.cfu_table != "NULL":
        cfu_table = pd.read_csv(args.cfu_table)
        cfu_dict = dict(zip(cfu_table.iloc[:, 0], cfu_table.iloc[:, 1]))
        print(f"Loaded CFU data for {len(cfu_dict)} samples")

    # Parse reference columns
    ref_cols = [int(x.strip()) for x in args.reference_columns.split(',')]
    print(f"Reference columns: {ref_cols}")

    # Calculate reference vector (sum across reference columns)
    reference_vec = freq_table.iloc[:, ref_cols].sum(axis=1).values
    print(f"Reference vector: {np.sum(reference_vec > 0)} non-zero barcodes")

    # Get sample columns (exclude reference)
    sample_cols = [i for i in range(len(freq_table.columns)) if i not in ref_cols]
    sample_names = freq_table.columns[sample_cols]
    print(f"Analyzing {len(sample_names)} samples: {list(sample_names)}\n")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Analyze each sample
    all_results = []

    for sample_name in sample_names:
        output_vec = freq_table[sample_name].values
        cfu = cfu_dict.get(sample_name, None)

        results = analyze_sample_clean(
            sample_name,
            reference_vec,
            output_vec,
            cfu=cfu,
            output_dir=output_dir,
            make_plots=not args.no_plots
        )

        if results is not None:
            all_results.append(results)

    # Create summary table
    if len(all_results) > 0:
        summary_df = pd.DataFrame([
            {
                'Sample': r['sample_name'],
                'Nb': r['Nb'],
                'Ns': r['Ns'],
                'F_stat': r['F_stat'],
                'Richness': r['richness'],
                'Shannon_H': r['shannon_diversity'],
                'Simpson': r['simpson_diversity'],
                'Evenness': r['evenness'],
                'AvgFrequency': r['avg_frequency'],
                'CFU': r['CFU'],
                'Log10_Ns': r['log10_Ns'],
                'Log10_CFU': r['log10_CFU'],
                'CFU_over_Ns': r['CFU_over_Ns'],
                'TotalCounts': r['total_counts']
            }
            for r in all_results
        ])

        # Save outputs
        summary_path = output_dir / 'FP_estimates_clean.csv'
        summary_df.to_csv(summary_path, index=False)

        print("\n" + "="*70)
        print("ANALYSIS COMPLETE!")
        print("="*70)
        print(f"\nResults saved to: {output_dir}")
        print(f"  - Summary table: {summary_path}")
        if not args.no_plots:
            print(f"  - Diagnostic plots: {output_dir}/*_FP_analysis.png")
        print("\n" + "="*70)
        print("\nSUMMARY TABLE:")
        print("="*70)

        # Print formatted summary
        display_cols = ['Sample', 'Nb', 'Ns', 'Richness', 'Shannon_H', 'CFU', 'CFU_over_Ns']
        print(summary_df[display_cols].to_string(index=False))
        print("="*70)
    else:
        print("\nNo results generated - check your input data!")


if __name__ == '__main__':
    main()
