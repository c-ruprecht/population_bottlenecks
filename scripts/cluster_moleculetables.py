### Clustering of UMIs in a molecule table using Gavage columns.
### Creates sums for all samples
from umi_tools import UMIClusterer


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process UMI clustering.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input moleculetables")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory")
    
    return parser.parse_args()


def cluster_umis(moleculetable, cluster_method = "directional", threshold = 1):
    df.set_index('umi_seq', inplace = True)
    #get gavage columns, create sum value and :
    df_gavage = df[[col for col in df.columns if 'gavage' in str(col).lower()]].copy()
    dict_gavage = df_gavage.sum(axis=1).to_dict()

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
    df_cluster = cluster_umis(df, cluster_method = "directional", threshold = 1)
    df_cluster.to_csv(args.output + '_clustered.csv')

    
if __name__ == "__main__":
    main()