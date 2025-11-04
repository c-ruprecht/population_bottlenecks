
# Create read tables to use with stampr
# Use a threshhold forkeeping a read using percentage cutoffs
# 1. Read total read tables and summary file to get total read counts
# conda run -n umi_tools_env python /Users/ruprec01/Documents/Faith_lab/Git/bc_seq_P4C2/scripts/create_molecule_table_R2.py -i /Volumes/sd/sample_data/processed/ -o /Volumes/sd/sample_data/moleculetables/ -g "s82" -t 500 -st "ST1,ST2,ST4,ST5,ST6"
import argparse
import os
import pandas as pd
import re
import plotly.express as px

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process UMI clustering.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input .fastq.gz file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory")
    parser.add_argument("-g", "--gavage_ids", required=True, help="Path to the gavage IDs file")
    parser.add_argument('-t', '--threshold', required=True, help='Threshold for normalization')
    parser.add_argument('-st', '--strain_barcode', required=False, help='Strain barcodes to keep')
    parser.add_argument('-l', '--locations', required==True, help='a tsv seperated sheet of samples and gavage locations')
    return parser.parse_args()

def get_readstable(path, path_out,li_gavage_ids, threshold, strain_barcode = None):
    li_gavage_ids = li_gavage_ids.split(',')

    # Initialize an empty DataFrame
    df_cage = pd.DataFrame()
    # Iterate over all subdirectories and files
    # only keep values over threshold
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.startswith('._'):
                continue
            if file.endswith('unique_umi_table.csv'):
                file_path = os.path.join(root, file)
                print(file_path)
                df = pd.read_csv(file_path)
                df['sample'] = file.split('_')[0]
                #make this dependent on gavage id
                if any(gavage_id in file for gavage_id in li_gavage_ids):
                    df_loc = df.loc[df['molecules']>1]
                else:
                    df_loc = df.loc[df['molecules']>int(threshold)]
                df_cage = pd.concat([df_cage, df_loc], ignore_index=True)
    
    print(df_cage)
    df_cage['strain'] = df_cage['index'].str.split('-').str[0]
    if strain_barcode is not None:
        df_cage = df_cage[df_cage['strain'].isin(strain_barcode.split(','))].copy()
    df_cage['molecules'] = df_cage['molecules'].astype(int)
    ## Total strain readtables, no mouse
    for strain in df_cage['strain'].unique():
        df_strain = df_cage[df_cage['strain'] == strain].copy()
        ###pivot table
        df_strain = df_strain.pivot_table(index='index', columns='sample', values='molecules', fill_value=0)
        print(df_strain)
        print(df_strain.columns)
        #resort to have inoculum first
        # sum all gavage columns
        # Add a column that sums all gavage columns
        gavage_columns = [col for col in df_strain.columns if any(gavage_id in col for gavage_id in li_gavage_ids)]
        if gavage_columns:  # Only create the sum column if there are gavage columns
            df_strain['s-gavage'] = df_strain[gavage_columns].sum(axis=1)
        sorted_columns =  ['s-gavage'] + \
                          [col for col in df_strain.columns if any(gavage_id in col for gavage_id in li_gavage_ids)] + \
                          [col for col in df_strain.columns if not any(gavage_id in col for gavage_id in li_gavage_ids) and col != 's-gavage']
        
        df_strain = df_strain[sorted_columns]
        # Convert all numeric columns to integers


        #fillna with 0s
        df_strain.fillna(0, inplace=True)
        numeric_columns = df_strain.select_dtypes(include=['number']).columns
        df_strain[numeric_columns] = df_strain[numeric_columns].astype(int)
        #export strain readstable
        df_strain.to_csv(path_out + '/' + strain + '_readstable.csv')

        #Create a figure to show gavage and overlay every other sample
        #df_strain = df_strain.reset_index()
        df_strain = df_strain.sort_values('s-gavage', ascending = True)
        for col in df_strain.columns:
            df_col = df_strain.loc[df_strain[col] > 0]
            fig = px.scatter(data_frame = df_col,
                            x = df_col.index,
                            y = ['s-gavage', col],
                            template = 'simple_white',
                            title = col)
            fig.write_image(path_out + '/figures/' + strain + '_'+col+'_scatter.png')


def main():
    args = parse_arguments()
    get_readstable(args.input, args.output, args.gavage_ids, args.threshold, args.strain_barcode)

if __name__ == "__main__":
    main()


