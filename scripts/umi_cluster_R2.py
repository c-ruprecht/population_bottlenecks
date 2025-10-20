# read fastq file
# test usage

#conda run -n umi_tools_env python /Users/ruprec01/Documents/Faith_lab/Git/bc_seq_P4C2/scripts/umi_cluster_R2.py -i /Volumes/sd/sample_data/trimmed/NTC-G5_R2_001.fastq.gz-trimmed.fastq -o /Volumes/sd/sample_data -s test -bc /Users/ruprec01/Documents/Faith_lab/Git/bc_seq_P4C2/input_files/strain_barcodes.fasta
import argparse
import gzip
from Bio import SeqIO
from umi_tools import UMIClusterer
import pandas as pd
import os
import plotly.express as px
import tempfile


# Set the temporary directory when running on cluster
os.makedirs("/sc/arion/scratch/ruprec01/tmp", exist_ok=True)
tempfile.tempdir = "/sc/arion/scratch/ruprec01/tmp"
os.environ["TMPDIR"] = "/sc/arion/scratch/ruprec01/tmp"

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process UMI clustering.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input .fastq.gz file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory")
    parser.add_argument("-s", "--sample", required=True, help="Sample name")
    parser.add_argument("-bc", "--strain_barcode", required=True, help="Path to the strain barcode fasta file")
    #parser.add_argument("-bc", "--strain_barcode", required=True, help="Path to the strain barcode fasta file")
    return parser.parse_args()

def process_fastq(file_path, len_barcode_region = 38, bc_dict_pre = {}, bc_dict_suff = {}):
    
    results = []
    # Handle gzipped files properly
    if file_path.endswith('.gz'):
        # Open gzipped file in text mode
        handle = gzip.open(file_path, 'rt')
    else:
        # Open regular file
        handle = open(file_path, 'r')
    #speed test version A:
    i = 0
    for record in SeqIO.parse(handle, "fastq"):
        sequence = str(record.seq[0:len_barcode_region])

        # remove all sequences that are not full length
        prefix = sequence[:10]
        suffix = sequence[-10:]
        umi_seq = sequence[10:-10]
        #maybe this shouuld be grabbed with 10:10+18

        # Map barcodes directly during parsing
        prefix_strain = bc_dict_pre.get(prefix)
        suffix_strain = bc_dict_suff.get(suffix)
        # Determine strain
        strain = None
        if prefix_strain == suffix_strain:
            strain = prefix_strain
        elif suffix_strain is None and prefix_strain is not None:
            strain = prefix_strain
        elif prefix_strain is None and suffix_strain is not None:
            strain = suffix_strain
            
        # Only keep reads with assigned barcodes and proper umi length
        if strain is not None and len(umi_seq) == 18:
            results.append({
                'id': record.id,
                'umi_seq': strain + '-' + umi_seq,
                'strain': strain,
            })
        i=i+1
    
    # Return as DataFrame directly
    return pd.DataFrame(results), i

def cluster_umis(umi_dict, cluster_method = "directional", threshold = 2):
    clusterer = UMIClusterer(cluster_method=cluster_method)
    
    #convert string keys to bytes
    umi_dict = {key.encode(): value for key, value in umi_dict.items()}
    clustered_umis = clusterer(umi_dict, threshold=threshold)
    li_clusters = [len(cluster) for cluster in clustered_umis]
    sum_clusters = sum(li_clusters)
    print('clustered '+ str(len(clustered_umis)) + ' of '+ str(sum_clusters) + ' total umis')
    # Create new dictionary with the clustered UMIs and total counts
    clustered_dict = {}
    #decode dict keys to str
    # Loop over the clustered UMIs
    for cluster in clustered_umis:
        for i, sub_cluster in enumerate(cluster):
            # get dict_umi value for each umi in cluster
            if i == 0:
                clustered_dict[cluster[0]] = +umi_dict[sub_cluster]

            else:
                clustered_dict[cluster[0]] = clustered_dict[cluster[0]]+umi_dict[sub_cluster]
    clustered_dict = {key.decode(): value for key, value in clustered_dict.items()}
    return clustered_dict



def main():
    args = parse_arguments()
    path = args.input 
    sample_name = args.sample

    # get barcode dictionaries
    bc_dict_pre = {}
    bc_dict_suff = {}
    #speed test version A:
    for record in SeqIO.parse(args.strain_barcode, "fasta"):
            umi = str(record.seq.reverse_complement())
            bc_dict_pre[umi[:10]] = record.id
            bc_dict_suff[umi[-10:]] = record.id
    print(bc_dict_pre)
    print(bc_dict_suff)
    #process fastq input files
    df_umis, tota_reads = process_fastq(args.input, len_barcode_region = 38, bc_dict_pre = bc_dict_pre, bc_dict_suff = bc_dict_suff)
    print(str(len(df_umis)) + " bardode reads of " + str(tota_reads) + ' reads processed')
    print(df_umis)
    #df_umis.to_csv(args.output + '/' + sample_name + '_umi_table.csv', index=False)
    df_gpd = df_umis.groupby(['strain', 'umi_seq']).size().reset_index(name='counts').sort_values(by='counts', ascending=False)
    df_gpd.to_csv(args.output + '/' + sample_name + '_unique_umi_table.csv', index=False)
    
    ### Run UMI clustering
    df_strains = pd.DataFrame()
    for strain in df_gpd['strain'].unique():
        df_loc= df_gpd.loc[df_gpd['strain'] == strain].copy()
        df_loc['idx'] = df_loc['umi_seq']
        df_loc = df_loc.set_index(['idx'])
        umi_dict = df_loc['counts'].to_dict()# Save umi_dict to the output file
        #adjacency, cluster, directionsal
        df_umis = pd.DataFrame.from_dict(cluster_umis(umi_dict, threshold=3,cluster_method='directional'), orient='index', columns=['count'])
        # reset index to get umi_seq as column
        df_umis = df_umis.reset_index()
        #df_umis['strain'] = file.split('-umi.fastq.gz')[0]
        df_umis = df_umis.sort_values(by='count', ascending=False)
        df_strains = pd.concat([df_strains, df_umis], axis=0)

    df_strains.to_csv(args.output + '/' + sample_name + '_total_readstable.csv', index=False)

    df_st3_spike = df_strains[df_strains['index'] == 'ST3-TCACACATTGACCTATGG'].copy()
    df_st3_spike['molecules'] = 1000
    df_st3_spike['reads/molecule'] = df_st3_spike['count'] / df_st3_spike['molecules']
    df_st3_spike

    df_strains['molecules'] = df_strains['count'] / df_st3_spike['reads/molecule'].values[0]
    df_strains['over_spike'] = df_strains['molecules'] >= 1000
    df_strains['strain'] = df_strains['index'].apply(lambda x: 'ST3-spike' if x == 'ST3-TCACACATTGACCTATGG' else x.split('-')[0] )
    df_strains['index'] = df_strains['index'].replace('ST3-TCACACATTGACCTATGG', 'ST3-spike')
    df_strains = df_strains.sort_values(by='molecules', ascending=False)
    
    #export molecules table
    df_strains[['index', 'molecules']].to_csv(args.output + '/' + sample_name + '_total_molecules_table.csv', index=False)

    df_strains = df_strains.sort_values(by=['strain', 'molecules'], ascending=True)
    # use defined color dictionary for STs
    color_dict = {'ST3-spike': 'black',
                    'ST3': 'grey',
                  'ST1': px.colors.qualitative.D3[0],
                  'ST2': px.colors.qualitative.D3[1],
                  'ST4': px.colors.qualitative.D3[2],
                  'ST5': px.colors.qualitative.D3[3],
                  'ST6': px.colors.qualitative.D3[4],
                  'ST7': px.colors.qualitative.D3[5],
                  'ST8': px.colors.qualitative.D3[6],
                  'ST9': px.colors.qualitative.D3[7]}
    
    fig = px.scatter(df_strains, x='index', 
                     y='molecules', 
                     color='strain', 
                     log_y=True,
                     color_discrete_map=color_dict,
                     width=2000, height=1000,
                    template='simple_white', 
                    title=args.sample,)
    fig.add_hline(y=1000, line_color="black", line_dash="dash", line_width=2)
    fig.write_image(args.output + '/' + sample_name + '_total_molecules_stable.png', width=2000, height=1000)    
if __name__ == "__main__":
    main()