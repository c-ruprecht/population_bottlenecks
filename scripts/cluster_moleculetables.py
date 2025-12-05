### Clustering of UMIs in a molecule table using Gavage columns.
### Filters 99 percentile of NTC or Water controls
### Creates clustered, denoised, and denoised without control files (sometimes NTC leads to getFP failures)
#/Volumes/sd/faith/MTCSB/projects/P4-barcoding_strains/20251104-analysis/test_empty/P4C1T8/test_clsuter
from umi_tools import UMIClusterer
import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process UMI clustering.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input moleculetables")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory")
    parser.add_argument("-t", "--threshold", required = False, help = " Set threshhold for umi cluster algorithm, default, directional 1")
    return parser.parse_args()


def cluster_umis(df, cluster_method = "directional", threshold = 1):
    df.set_index('umi_seq', inplace = True)
    #get gavage columns, create sum value
    df_gavage = df[[col for col in df.columns if 'gavage' in str(col).lower()]].copy()
    dict_gavage = df_gavage.sum(axis=1).to_dict()
    dict_gavage = {k: v for k, v in dict_gavage.items() if v > 0}

    # Encode keys for umi_tools cluster
    umi_dict = {key.encode(): value for key, value in dict_gavage.items()}
    clusterer = UMIClusterer(cluster_method=cluster_method)
    clustered_umis = clusterer(umi_dict, threshold=threshold)
    li_clusters = [len(cluster) for cluster in clustered_umis]
    sum_clusters = sum(li_clusters)
    print('clustered '+ str(len(clustered_umis)) + ' of '+ str(sum_clusters) + ' total umis')

    # Build a list of Series, then concat once at the end
    cluster_series_list = []

    for cluster in clustered_umis:
        # Get representative UMI (first in cluster)
        representative_umi = cluster[0].decode()
        
        # Get all UMIs in this cluster
        cluster_umis = [umi.decode() for umi in cluster]
        
        # Get all rows from original dataframe for UMIs in this cluster
        cluster_rows = df.loc[cluster_umis]
        
        # Sum ALL columns across these rows
        summed_row = cluster_rows.sum(axis=0)
        summed_row.name = representative_umi  # Set the index name
        
        # Append to list instead of adding to DataFrame
        cluster_series_list.append(summed_row)

    # Concat all at once - much faster!
    clustered_df = pd.concat(cluster_series_list, axis=1).T
    clustered_df.reset_index(inplace = True)
    clustered_df = clustered_df.rename(columns = {'index': 'umi_seq'})
    return clustered_df


def main():
    args = parse_arguments()

    df = pd.read_csv(args.input)
    print('starting to cluster')
    print(df)
    df_cluster = cluster_umis(df, cluster_method = "directional", threshold = args.threshold)

    #downscale reads if neccessary to stay in integer R limit
    SAFE_THRESHOLD = 2e9
    numeric_cols = df_cluster.select_dtypes(include=['int64', 'float64']).columns
    max_sample_total = df_cluster[numeric_cols].sum(axis=0).max()
    print(f"  Max molecules per sample: {max_sample_total:,.0f}")
    if max_sample_total > SAFE_THRESHOLD:
        scale_factor = max_sample_total / SAFE_THRESHOLD
        print(f"  WARNING: Exceeds safe threshold!")
        print(f"  Scaling down by factor: {scale_factor:.2f}")
        
        # Downsample all numeric columns proportionally
        df_cluster[numeric_cols] = (df_cluster[numeric_cols] / scale_factor).round().astype(int)
        
        # Verify new max
        new_max = df_cluster[numeric_cols].sum(axis=0).max()
        print(f"  New max sample total: {new_max:,.0f}")
    
    # Convert numeric columns to integers to save space
    df_cluster[numeric_cols] = df_cluster[numeric_cols].astype(int)
    
    #drop everything where gavage.lower() is 0
    gavage_cols = [col for col in df_cluster.columns if 'gavage' in col.lower()]
    df_cluster_filtered = df_cluster[df_cluster[gavage_cols].sum(axis=1) > 0].copy()
    df_cluster_filtered.to_csv(args.output + '_clustered.csv', index = False)
    print('exported clustered file')

    ### Noise filter based on negative controls 99 percentile
    #ntc_cols = [col for col in df_cluster_filtered.columns if 'NTC' in col or 'WATER' in col.upper()]

    #if len(ntc_cols) > 0:
    #    df_stack = df_cluster_filtered.set_index('umi_seq')[ntc_cols].stack().reset_index()
    ##    df_stack.columns = ['umi_seq', 'sample', 'clustered_molecules']
    #    noise_threshold = df_stack['clustered_molecules'].quantile(0.99)
    #    print(f'Noise threshold (99th percentile): {noise_threshold}')

    #    df_denoised = df_cluster_filtered.copy()
    #    numeric_cols = df_denoised.select_dtypes(include=['int64', 'float64']).columns
    #    df_denoised[numeric_cols] = df_denoised[numeric_cols].where(df_denoised[numeric_cols] >= noise_threshold, 0)
   # else:
    #    print('No NTC or WATER columns found, skipping denoising')
   #     df_denoised = df_cluster_filtered.copy()

    #df_denoised.to_csv(args.output + '_clustered_denoised.csv', index = False)
    #print('exported denoised file')

    # Create version without control columns
    ##cols_to_keep = ['umi_seq'] + [col for col in df_denoised.columns
    #                                if col != 'umi_seq'
    #                                and 'NTC' not in col
    #                                and 'WATER' not in col.upper()]
    #df_denoised_no_controls = df_denoised[cols_to_keep].copy()
    #df_denoised_no_controls.to_csv(args.output + '_clustered_denoised_nocontrols.csv', index = False)
    #print('exported denoised without controls file')

    #print('All output files created successfully')
    return 0

if __name__ == "__main__":
    import sys
    try:
        sys.exit(main())
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)