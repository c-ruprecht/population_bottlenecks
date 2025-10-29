"""
getFP - Founder Population Analysis for Barcode Sequencing Data

This script estimates the founding population size (Ns) and bottleneck population (Nb)
from barcode sequencing experiments using frequency-based statistical methods.

Usage:
    python getFP.py <reads_table.csv> <cfu_table.csv|NULL> <reference_columns> <min_weight> <output_dir>

Example:
    python getFP.py data/readstable.csv NULL "0" 0.5 output/
"""

import numpy as np
import pandas as pd
import argparse
from scipy import stats
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
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
    ratio = np.sum(output_vec > 1) / max(np.sum(output_vec == 1), 1)

    # If singletons dominate (ratio < 1.2), apply noise correction
    if np.sum(output_vec == 1) > 1.2 * np.sum(output_vec > 1):
        output_vec_corrected = output_vec - 1
        output_vec_corrected[output_vec_corrected < 0] = 0
        return output_vec_corrected

    return output_vec


def find_greatest_frequency_drop(output_vec):
    """
    Find the largest drop in log-frequency space.

    Large drops indicate transitions between real populations and noise.
    """
    output_copy = output_vec.copy()
    output_copy[output_copy < 2] = 0

    sorted_counts = np.sort(output_copy[output_copy > 0])[::-1]

    if len(sorted_counts) < 2:
        return None

    # Calculate log-differences
    shifted = np.append(sorted_counts[1:], sorted_counts[-1])
    log_diff = np.log(sorted_counts + 1) - np.log(shifted + 1)
    log_diff = log_diff[:-1]

    # Threshold: log(10) ≈ 2.3 means 10-fold drop
    if len(log_diff) > 0 and np.max(log_diff) > 2.302585:
        return np.argmax(log_diff)

    return None


def calculate_resiliency_curve(input_vec, output_vec, max_iterations=20000):
    """
    Calculate the "resiliency curve" - bottleneck size estimates across barcode ranks.

    This iteratively estimates Nb by progressively removing the least abundant barcodes,
    revealing transitions between true populations and noise.

    Parameters:
    -----------
    input_vec : array
        Reference population counts
    output_vec : array
        Sample counts
    max_iterations : int
        Maximum number of barcodes to analyze

    Returns:
    --------
    resiliency_values : array
        Nb estimates at each iteration
    """
    # Sort by abundance (least to most)
    sorted_indices = np.argsort(output_vec)

    input_sorted = input_vec[sorted_indices]
    output_sorted = output_vec[sorted_indices]

    n_barcodes = len(input_vec)
    iterations = min(max_iterations, n_barcodes)

    resiliency_values = []

    # Calculate Nb for each subset (removing least abundant first)
    for i in range(n_barcodes - iterations, n_barcodes):
        subset_input = input_sorted[:i+1]
        subset_output = output_sorted[:i+1]

        if np.sum(subset_output) == 0 or np.sum(subset_input) == 0:
            resiliency_values.append(1)
            continue

        Nb, _ = calculate_wright_fisher_bottleneck(subset_input, subset_output)
        resiliency_values.append(Nb)

    return np.array(resiliency_values)


def find_population_breaks_mcmc(resiliency_curve, n_samples=10000):
    """
    Find breaks in resiliency curve using Monte Carlo random walk.

    This uses a stochastic search to find local minima in the resiliency curve,
    which represent transitions between distinct barcode populations.

    Parameters:
    -----------
    resiliency_curve : array
        Nb estimates across iterations
    n_samples : int
        Number of MCMC samples

    Returns:
    --------
    break_points : array
        Indices of detected population breaks
    """
    if len(resiliency_curve) == 0 or np.nanmean(resiliency_curve) == 1:
        return [len(resiliency_curve)]

    # Sample starting points logarithmically
    n_points = len(resiliency_curve)
    sample_indices = np.linspace(0, n_points - 1,
                                  num=int(np.log(n_points) + 1),
                                  dtype=int)

    break_points = []

    for start_idx in sample_indices:
        if start_idx >= len(resiliency_curve):
            continue

        current_value = resiliency_curve[start_idx]
        current_min = current_value

        # Random walk to find local minimum
        for _ in range(n_samples):
            # Generate new position with normal distribution
            sd = max(np.sum(resiliency_curve > 1) / 20, 1)
            new_idx = int(abs(np.random.normal(start_idx, sd)))

            # Boundary checks
            if new_idx >= len(resiliency_curve) or new_idx < 0:
                continue

            new_value = resiliency_curve[new_idx]
            if new_value < current_min:
                current_min = new_value
                start_idx = new_idx

        # Find first occurrence of minimum
        min_indices = np.where(resiliency_curve == current_min)[0]
        if len(min_indices) > 0:
            break_points.append(min_indices[0])

    return np.unique(break_points)


def calculate_cumulative_fraction(output_vec, n_top_barcodes):
    """Calculate what fraction of reads the top N barcodes account for."""
    sorted_counts = np.sort(output_vec)[::-1]
    top_counts = sorted_counts[:n_top_barcodes]
    return np.sum(top_counts) / np.sum(output_vec)


def estimate_founding_population_simulation(input_vec, output_vec, n_detected_barcodes):
    """
    Estimate founding population size (Ns) using rarefaction simulation.

    This simulates how many barcodes would be detected at different founding
    population sizes and interpolates to find Ns that matches observations.

    Parameters:
    -----------
    input_vec : array
        Reference population
    output_vec : array
        Observed sample
    n_detected_barcodes : int
        Number of detected barcodes (after noise filtering)

    Returns:
    --------
    Ns : float
        Estimated founding population size
    interpolation_curve : tuple
        (x_values, y_values) for the rarefaction curve
    """
    # First resample: create expected distribution based on reference
    total_reads = np.sum(output_vec)
    reference_freq = input_vec / np.sum(input_vec)

    # Multinomial resample to get expected counts
    first_resample = np.random.multinomial(total_reads, reference_freq)

    # Rarefaction: how many unique barcodes at different founding sizes?
    n_unique_barcodes = np.sum(input_vec != 0)

    # Sample at different depths
    steps_low = np.linspace(1, n_unique_barcodes, 100, dtype=int)
    steps_high = np.linspace(n_unique_barcodes, n_unique_barcodes * 20, 100, dtype=int)
    steps = np.concatenate([steps_low, steps_high])

    expected_detections = []

    for n_founders in steps:
        # Simulate sampling: which barcodes would be detected?
        n_simulations = 5
        detections = []

        for _ in range(n_simulations):
            sampled = np.random.multinomial(n_founders,
                                           first_resample / np.sum(first_resample))
            n_detected = np.sum(sampled > 0)
            detections.append(n_detected)

        expected_detections.append(np.mean(detections))

    # Interpolate to find Ns
    interp_func = interp1d(steps, expected_detections,
                           kind='linear',
                           fill_value='extrapolate')

    # Find founding size that gives observed number of barcodes
    x_fine = np.linspace(np.min(steps), np.max(steps), len(input_vec))
    y_fine = interp_func(x_fine)

    # Find closest match
    closest_idx = np.argmin(np.abs(y_fine - n_detected_barcodes))
    Ns = x_fine[closest_idx]

    return Ns, (x_fine, y_fine)


def determine_noise_threshold(indices_df, min_weight=0.5):
    """
    Determine where noise begins in the barcode distribution.

    Uses weight-based thresholding to separate true populations from noise.

    Parameters:
    -----------
    indices_df : DataFrame
        Columns: ['n_barcodes', 'fraction_reads']
    min_weight : float
        Minimum fraction of reads for a population to be considered real

    Returns:
    --------
    noise_start : int
        Number of barcodes before noise begins
    """
    # Find largest weight drop
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

    return noise_start


def analyze_sample(sample_name, input_vec, output_vec, cfu=None,
                   min_weight=0.5, output_dir='.', make_plots=True):
    """
    Main analysis function for a single sample.

    Performs complete FP (Founder Population) analysis including:
    1. Noise correction
    2. Population structure detection
    3. Bottleneck size (Nb) estimation
    4. Founding population (Ns) estimation

    Parameters:
    -----------
    sample_name : str
        Sample identifier
    input_vec : array
        Reference population barcode counts
    output_vec : array
        Sample barcode counts
    cfu : float, optional
        Colony forming units (if available)
    min_weight : float
        Minimum population weight threshold
    output_dir : str
        Output directory path
    make_plots : bool
        Whether to generate diagnostic plots

    Returns:
    --------
    results : dict
        Dictionary containing all estimated parameters
    """
    print(f"\nAnalyzing sample: {sample_name}")

    # Initial cleanup
    total_reads = np.sum(output_vec)
    output_vec = remove_low_frequency_noise(output_vec)

    # Handle CFU
    if cfu is None or cfu == 0 or (isinstance(cfu, float) and np.isnan(cfu)):
        cfu = 1e20

    # Determine iteration limit based on singleton ratio
    ratio = np.sum(output_vec > 1) / max(np.sum(output_vec == 1), 1)
    factor = max(0, 1.5 - ratio)
    remove_ones = min(np.sum(output_vec == 1),
                      round(np.sum(output_vec == 1) * factor))

    if np.sum(output_vec > 1) > 1.5 * np.sum(output_vec == 1):
        max_iterations = min(int(np.ceil(cfu)), np.sum(output_vec > 0))
    else:
        max_iterations = min(int(np.ceil(cfu)),
                            np.sum(output_vec > 0) - remove_ones)

    # Noise correction
    output_vec = adjust_singleton_noise(output_vec)
    greatest_freq_drop = find_greatest_frequency_drop(output_vec)

    # Calculate resiliency curve
    resiliency_curve = calculate_resiliency_curve(input_vec, output_vec, max_iterations)

    # Find population breaks
    break_points = find_population_breaks_mcmc(resiliency_curve)

    # Add frequency-based breaks
    if greatest_freq_drop is not None:
        break_points = np.append(break_points, greatest_freq_drop)

    # Add resiliency curve breaks
    if len(resiliency_curve) > 1:
        log_res = np.log(resiliency_curve + 1)
        res_diff = log_res[1:] - log_res[:-1]
        if len(res_diff) > 0 and np.max(res_diff) > 2.302585:
            greatest_res_break = np.argmax(res_diff)
            break_points = np.append(break_points, greatest_res_break)

    # Finalize breaks
    break_points = np.append(break_points, len(resiliency_curve))
    break_points = np.unique(break_points[break_points > 0])
    break_points = np.sort(break_points)

    # Remove adjacent breaks
    if len(break_points) > 1:
        diff_breaks = np.diff(break_points)
        if np.min(diff_breaks) == 1 and cfu > 2:
            to_remove = np.where(diff_breaks == 1)[0]
            if len(to_remove) > 0:
                break_points = np.delete(break_points, to_remove[-1])

    # Calculate population fractions
    fractions = []
    for n_barcodes in break_points:
        frac = calculate_cumulative_fraction(output_vec, n_barcodes)
        fractions.append(frac)

    # Create indices table
    staggered = np.concatenate([[0], fractions[:-1]])
    subtracted = np.array(fractions) - staggered
    subtracted = subtracted / np.sum(subtracted)

    indices_df = pd.DataFrame({
        'n_barcodes': break_points,
        'fraction_reads': subtracted
    })

    # Determine noise threshold
    noise_start = determine_noise_threshold(indices_df, min_weight)

    # Filter noise
    output_vec_clean = output_vec.copy()
    n_zeros = len(input_vec) - noise_start
    sorted_indices = np.argsort(output_vec_clean)
    output_vec_clean[sorted_indices[:n_zeros]] = 0

    # Calculate bottleneck size
    Nb, F_stat = calculate_wright_fisher_bottleneck(input_vec, output_vec)

    # Calculate founding population size
    Ns, interp_curve = estimate_founding_population_simulation(
        input_vec, output_vec, noise_start
    )

    # Average frequency
    avg_freq = 1 / geometric_mean(
        output_vec_clean[output_vec_clean > 0] / np.sum(output_vec_clean)
    )

    # Ns at minweight cutoff
    reads_at_cutoff = (1 - min_weight) * np.sum(output_vec)
    sorted_reads = np.sort(output_vec[output_vec > 0])
    cumsum_reads = np.cumsum(sorted_reads)

    if len(cumsum_reads[cumsum_reads > reads_at_cutoff]) > 0:
        min_cutoff_idx = np.where(cumsum_reads > reads_at_cutoff)[0][0]
        min_cutoff = sorted_reads[min_cutoff_idx]
        n_barcodes_at_minweight = np.sum(output_vec > min_cutoff) + 1

        # Interpolate Ns at this cutoff
        x_vals, y_vals = interp_curve
        closest = np.argmin(np.abs(y_vals - n_barcodes_at_minweight))
        Ns_min_cutoff = x_vals[closest]
    else:
        Ns_min_cutoff = Ns

    # Plotting
    if make_plots:
        create_diagnostic_plots(
            sample_name, output_vec, resiliency_curve,
            break_points, indices_df, Nb, Ns, avg_freq,
            Ns_min_cutoff, output_dir
        )

    # Prepare results
    if cfu == 1e20:
        cfu = 0

    results = {
        'sample_name': sample_name,
        'total_reads': total_reads,
        'n_barcodes': noise_start,
        'Ns_min_cutoff': Ns_min_cutoff,
        'Nb': Nb,
        'Ns': Ns,
        'avg_frequency': avg_freq,
        'CFU': cfu,
        'log10_Ns': np.log10(Ns) if Ns > 0 else 0,
        'log10_CFU': np.log10(cfu) if cfu > 0 else 0,
        'CFU_over_Ns': cfu / Ns if Ns > 0 else 0,
        'filtered_counts': output_vec_clean
    }

    print(f"Results for {sample_name}:")
    print(f"  Nb = {Nb:.2f}")
    print(f"  Ns = {Ns:.2f}")
    print(f"  Average Frequency = {avg_freq:.2f}")
    print(f"  Ns at minweight = {Ns_min_cutoff:.2f}")

    return results


def create_diagnostic_plots(sample_name, output_vec, resiliency_curve,
                           break_points, indices_df, Nb, Ns, avg_freq,
                           Ns_min_cutoff, output_dir):
    """Create diagnostic plots showing analysis results."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))

    # Plot 1: Resiliency curve
    ax1 = axes[0]
    ax1.plot(range(1, len(resiliency_curve) + 1), resiliency_curve)
    ax1.set_yscale('log')
    ax1.set_ylim(0.01, 2e6)
    ax1.set_title('Resiliency Curve')
    ax1.set_ylabel('Nb')
    ax1.set_xlabel('Iteration')

    # Add break point lines
    colors = plt.cm.rainbow(np.linspace(0, 1, len(break_points)))
    for i, bp in enumerate(break_points):
        ax1.axvline(bp, color=colors[i], linewidth=1.5, alpha=0.7)
        ax1.text(bp, np.max(resiliency_curve) * 0.8,
                f'G{i+1}: {bp}', rotation=90,
                verticalalignment='bottom', color=colors[i], fontsize=8)

    # Plot 2: Frequency distribution
    ax2 = axes[1]
    frequencies = output_vec / np.sum(output_vec)
    ax2.plot(frequencies)
    ax2.set_yscale('log')
    ax2.set_ylim(1e-8, 1)
    ax2.set_title(f'{sample_name} - Barcode Frequencies')
    ax2.set_ylabel('Frequency')
    ax2.set_xlabel('Barcode')

    # Add cutoff lines
    sorted_output = np.sort(output_vec)[::-1]
    for i, bp in enumerate(break_points):
        if bp < len(sorted_output):
            cutoff_freq = sorted_output[bp] / np.sum(output_vec)
            ax2.axhline(cutoff_freq, color=colors[i],
                       linewidth=1.5, alpha=0.7)
            ax2.text(len(output_vec) * 0.9, cutoff_freq * 1.2,
                    f'C{i+1}: {cutoff_freq:.2e}',
                    color=colors[i], fontsize=8)

    # Plot 3: Statistics text
    ax3 = axes[2]
    ax3.axis('off')
    stats_text = f"""
    Statistics:
    Nb = {Nb:.2f}
    Ns = {Ns:.2f}
    Average Frequency = {avg_freq:.2f}
    Ns at minweight = {Ns_min_cutoff:.2f}

    Population Breakpoints:
    {indices_df.to_string(index=False)}
    """
    ax3.text(0.1, 0.5, stats_text, fontsize=10,
            verticalalignment='center', family='monospace')

    plt.tight_layout()
    output_path = Path(output_dir) / f'{sample_name}_plot.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved plot to {output_path}")


def main():
    """Main entry point for command-line execution."""
    parser = argparse.ArgumentParser(
        description='FP Analysis: Estimate founding population from barcode sequencing'
    )
    parser.add_argument('reads_table', type=str,
                       help='Path to reads table CSV (barcodes × samples)')
    parser.add_argument('cfu_table', type=str,
                       help='Path to CFU table CSV or "NULL"')
    parser.add_argument('reference_columns', type=str,
                       help='Comma-separated column indices for reference (0-indexed)')
    parser.add_argument('min_weight', type=float,
                       help='Minimum weight threshold (e.g., 0.5)')
    parser.add_argument('output_dir', type=str,
                       help='Output directory path')

    args = parser.parse_args()

    # Load data
    reads_table = pd.read_csv(args.reads_table, index_col=0)

    # Load CFU table if provided
    cfu_dict = {}
    if args.cfu_table != "NULL":
        cfu_table = pd.read_csv(args.cfu_table)
        # Convert CFU column to numeric, coercing errors (like "-") to NaN
        cfu_values = pd.to_numeric(cfu_table.iloc[:, 1], errors='coerce')
        cfu_dict = dict(zip(cfu_table.iloc[:, 0], cfu_values))

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
    filtered_table = pd.DataFrame(index=reads_table.index)

    for sample_name in sample_names:
        output_vec = reads_table[sample_name].values
        cfu = cfu_dict.get(sample_name, None)

        results = analyze_sample(
            sample_name, reference_vec, output_vec,
            cfu=cfu, min_weight=args.min_weight,
            output_dir=output_dir, make_plots=True
        )

        all_results.append(results)
        filtered_table[sample_name] = results['filtered_counts']

    # Create summary table
    summary_df = pd.DataFrame([
        {
            'sample': r['sample_name'],
            'TotalReads': r['total_reads'],
            'Number of barcodes': r['n_barcodes'],
            'Ns_MinCutoff': r['Ns_min_cutoff'],
            'Nb': r['Nb'],
            'Ns': r['Ns'],
            'AverageFrequency': r['avg_frequency'],
            'CFU': r['CFU'],
            'Log10Ns': r['log10_Ns'],
            'Log10CFU': r['log10_CFU'],
            'CFU/Ns': r['CFU_over_Ns']
        }
        for r in all_results
    ])

    # Save outputs
    summary_df.to_csv(output_dir / 'TableOfEstimates.csv', index=False)
    filtered_table.to_csv(output_dir / 'FrequenciesWithoutNoise.csv')

    print("\n" + "="*60)
    print("Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print("="*60)
    print("\nSummary:")
    print(summary_df.to_string(index=False))


if __name__ == '__main__':
    main()
