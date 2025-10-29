
import pandas as pd
import numpy as np
import argparse
import os
from typing import Dict, List, Tuple

def get_genetic_distance(f1: np.ndarray, f2: np.ndarray) -> float:
    """Calculate genetic distance between two frequency vectors using chord distance."""
    # Normalize vectors
    f1 = f1 / np.sum(f1) if np.sum(f1) > 0 else f1
    f2 = f2 / np.sum(f2) if np.sum(f2) > 0 else f2
    
    # Calculate cosine similarity
    cos_theta = np.sum(np.sqrt(f1 * f2))
    
    # Correct for potential numerical errors
    if 1 - cos_theta < 0:
        cos_theta = 1
    
    # Calculate chord distance
    chord_distance = 2 * np.sqrt(2) / np.pi * np.sqrt(1 - cos_theta)
    
    return chord_distance

def calculate_resilience_distance(vec1: np.ndarray, vec2: np.ndarray, 
                                vec1_name: str, vec2_name: str, limit: int) -> Tuple:
    """Calculate resilience distance between two samples."""
    print(f"{vec1_name} {vec2_name}")
    
    # Calculate number of iterations
    times = min(np.sum(vec1 > 0), np.sum(vec2 > 0), limit)
    
    # Calculate geometric mean of frequencies
    vec1_freq = vec1 / np.sum(vec1) if np.sum(vec1) > 0 else vec1
    vec2_freq = vec2 / np.sum(vec2) if np.sum(vec2) > 0 else vec2
    square_root = np.sqrt(vec1_freq * vec2_freq)
    
    # Create a dataframe for sorting
    bind_df = pd.DataFrame({
        'vec1': vec1, 
        'vec2': vec2, 
        'square_root': square_root
    })
    
    # Sort by geometric mean (ascending)
    bind_sorted = bind_df.sort_values('square_root')
    
    # Iteratively remove strains and calculate genetic distance
    g_values = [0]  # Start with 0 as in R script
    
    for _ in range(times):
        current_vec1 = bind_sorted['vec1'].values
        current_vec2 = bind_sorted['vec2'].values
        gd = get_genetic_distance(current_vec1, current_vec2)
        g_values.append(gd)
        
        # Remove the largest strain
        if len(bind_sorted) > 1:
            bind_sorted = bind_sorted.iloc[:-1]
    
    # Remove the first element (initialization)
    g_values = g_values[1:]
    
    # Calculate RD (resilience distance) and initial GD
    rd = np.sum(np.array(g_values) < 0.8)
    gd = g_values[0] if g_values else 0
    
    return rd, gd, vec1_name, vec2_name

def correct_rd_values(rd_matrix: pd.DataFrame, reads_table: pd.DataFrame, limit: int) -> pd.DataFrame:
    """Correct RD values for limited iterations."""
    sample_names = rd_matrix.columns
    corrected_rd = pd.DataFrame(index=rd_matrix.index, columns=rd_matrix.columns)
    
    for name in sample_names:
        rd_vec = rd_matrix[name].values
        positions = np.where(rd_vec == limit)[0]
        
        if len(positions) == 0:
            corrected_rd[name] = rd_vec
        else:
            # Count non-zero values in sample
            n_barcodes = np.sum(reads_table[name].values != 0)
            rd_vec_corrected = rd_vec.copy()
            rd_vec_corrected[rd_vec == limit] = n_barcodes
            corrected_rd[name] = rd_vec_corrected
    
    return corrected_rd

def convert_to_frd(vec: np.ndarray) -> np.ndarray:
    """Convert RD values to FRD (Fractional Resilience Distance)."""
    vec = np.array(vec)
    vec = vec[~np.isnan(vec)]
    return np.log(vec + 1) / np.log(np.max(vec) + 1)

def resiliency_genetic_distance(reads_table_name: str, exp_name: str, 
                               comparison_metadata_name: str, limit: int, 
                               output_dir: str) -> Dict:
    """Main function to calculate resiliency genetic distance metrics."""
    # Read input files
    reads_table = pd.read_csv(reads_table_name, index_col=0)
    metadata = pd.read_csv(comparison_metadata_name)
    # Print original column names for debugging
    print("Original reads table columns:", list(reads_table.columns))
    print("Original metadata sample names:", list(metadata['Sample']))
    metadata['Sample'] = metadata['Sample'].str.replace('.', '-')  # THIS LINE IS MISSING
    reads_table.columns = [col.replace('.', '-') for col in reads_table.columns]
    # Validate metadata
    if "" in metadata['Sample'].values:
        print("STOP!!!")
        return {}
    
    # Get ordered sample names
    metadata = metadata.sort_values(['Group', 'Order'])
    sample_names = metadata['Sample'].values
    
    # Create all pairwise comparisons
    vector1 = np.repeat(sample_names, len(sample_names))
    vector2 = np.tile(sample_names, len(sample_names))
    table_of_comparisons = pd.DataFrame({'vector1': vector1, 'vector2': vector2})
    
    # Filter to comparisons within same group
    group_table = pd.DataFrame(index=table_of_comparisons.index, columns=['group1', 'group2'])
    
    for idx, row in table_of_comparisons.iterrows():
        group1 = metadata[metadata['Sample'] == row['vector1']]['Group'].values[0]
        group2 = metadata[metadata['Sample'] == row['vector2']]['Group'].values[0]
        group_table.loc[idx] = [group1, group2]
    
    filtered_comparisons = table_of_comparisons[group_table['group1'] == group_table['group2']]
    
    # Calculate RD for each comparison
    rd_results = []
    for _, row in filtered_comparisons.iterrows():
        vec1_name, vec2_name = row['vector1'], row['vector2']
        vec1 = reads_table[vec1_name].values
        vec2 = reads_table[vec2_name].values
        
        rd, gd, _, _ = calculate_resilience_distance(vec1, vec2, vec1_name, vec2_name, limit)
        rd_results.append([rd, gd, vec1_name, vec2_name])
    
    rd_table = pd.DataFrame(rd_results, columns=['RD', 'GD', 'vector1', 'vector2'])
    
    # Initialize matrices
    rd_matrix = pd.DataFrame(0, index=sample_names, columns=sample_names)
    gd_matrix = pd.DataFrame(0, index=sample_names, columns=sample_names)
    
    # Fill matrices
    for _, row in rd_table.iterrows():
        rd_matrix.loc[row['vector1'], row['vector2']] = row['RD']
        gd_matrix.loc[row['vector1'], row['vector2']] = row['GD']
    
    # Correct RD values
    corrected_rd_matrix = correct_rd_values(rd_matrix, reads_table, limit)
    
    # Calculate FRD matrix
    corrected_frd_matrix = pd.DataFrame(index=sample_names, columns=sample_names)
    for col in corrected_rd_matrix.columns:
        corrected_frd_matrix[col] = convert_to_frd(corrected_rd_matrix[col])
    
    # Write output files
    os.makedirs(output_dir, exist_ok=True)
    corrected_rd_matrix.to_csv(os.path.join(output_dir, f"CorrectedRD_{exp_name}.csv"))
    gd_matrix.to_csv(os.path.join(output_dir, f"GD_{exp_name}.csv"))
    corrected_frd_matrix.to_csv(os.path.join(output_dir, f"CorrectedFRD_{exp_name}.csv"))
    
    return {
        'RD_table': rd_table,
        'RD_matrix': rd_matrix,
        'Corrected_RD_matrix': corrected_rd_matrix,
        'GD_matrix': gd_matrix,
        'Corrected_FRD_matrix': corrected_frd_matrix
    }



def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Calculate resiliency genetic distance.")
    parser.add_argument("-rt", "--readstable", required=True, help="Path to reads table file")
    parser.add_argument("-n", "--name", required=True, help="Experiment name")
    parser.add_argument("-m", "--meta", required=True, help="Path to metadata table")
    parser.add_argument("-l", "--limit", required=True, type=int, help="Limit for iterations")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    return parser.parse_args()

def main() -> None:
    """Main function to run the script."""
    args = parse_arguments()
    
    results = resiliency_genetic_distance(
        args.readstable,
        args.name,
        args.meta,
        args.limit,
        args.output
    )
    
    print("Processing complete. Results saved to:", args.output)

if __name__ == "__main__":
    main()