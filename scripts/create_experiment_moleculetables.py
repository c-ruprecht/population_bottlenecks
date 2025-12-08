import argparse
import pandas as pd
from pathlib import Path



def parse_arguments():
    parser = argparse.ArgumentParser(description="Process UMI clustering.")
    parser.add_argument("-b", "--basepath", required = True, help = "Path to the ILL readtable folders")
    parser.add_argument("-m", "--metadir", required=True, help="Path to the metadata directories")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory")
    parser.add_argument('-e', "--empty", required=False, help="Path to a csv sheet containing characterizations of samples as empty")
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    ## read empty
    df_empty = pd.read_csv(args.empty)
    df_empty = df_empty[['ILL', 'sample', 'is_empty']].copy()
    df_empty.columns = ['ill', 'file_name','is_empty']
    df_empty['ill'] = 'ILL' + df_empty['ill'].astype(str)
    ### P4C3T1
    print(50* '=')
    print(' Creating moleculetables for P4C3T2 100K barcodese')
    print(50*'-')
    print(' ')
    try:
        basepath = Path(args.basepath)
        exportpath = Path(args.output) / 'P4C3T1/moleculetables'
        exportpath.mkdir(parents=True, exist_ok=True)

        ## Read metadata & ill sheets
        metadir = Path(args.metadir)
        #P4C2T4:
        df_ill = pd.read_excel(metadir / 'ill_sheets' / 'ILL134_sample_sheet.xlsx',
                                    sheet_name = 'Individual Library Information', header = 1)
        df_samplesheet = pd.read_csv(metadir / 'metadata' / 'P4C3T1_cr_mouse_sample - P4C3T1.tsv', sep = '\t')

        df_samplesheet = df_samplesheet.loc[df_samplesheet['sample_type'].str.contains('DNA')].copy()
        df_ill = df_ill.dropna(subset=['Sample Name/Pool Name*'])
        df_ill['sample_id']= df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).split('_')[0])
        df_ill['project'] = df_ill['Library Name*'].apply(lambda x: str(x).split('_')[1])
        #correct library name in ILL138 run
        df_ill['ill'] = df_ill['Library Name*'].apply(lambda x:  str(x).split('_')[0])
        df_ill['file_name'] = df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).replace('_', '-'))
        df_merge = pd.merge(df_ill[['sample_id', 'ill', 'file_name', 'project']], df_samplesheet, how = 'right', on = ['sample_id', 'project'])
        df_merge.head(10)

        df_ill_ntc = df_ill.loc[df_ill['Sample Name/Pool Name*'].str.contains("WATER|NTC")]
        df_merge = df_merge.dropna(subset= ['file_name'])
        print(df_merge)
        print(df_ill_ntc)

        
        li_strains = ['ST1', 'ST2']
        df_ntc = pd.DataFrame()
        for idx, row in df_ill_ntc.iterrows():
            file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
            df = pd.read_csv(file_path)
            df['sample'] = row['Sample Name/Pool Name*']+'_'+row['sample_id']
            df_ntc = pd.concat([df_ntc, df])
        
        df_ntc_pivot = df_ntc.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
        print('FOUND GAVAGE DOSES:')
        print(df_merge['gavage_dose_200ul'].unique())
        for dose in df_merge['gavage_dose_200ul'].unique():
            print(dose)
            df_dose = df_merge.loc[df_merge['gavage_dose_200ul']==dose]
            #get gavage reads dataframe
            #drop files without sequencing data
            df_merge.dropna(subset=['file_name'], inplace = True)
            df_gavage = df_dose[df_dose['sample_type'] == 'Gavage DNA'].copy()
            print(df_gavage)
            df_reads_gavage = pd.DataFrame()
            for idx, row in df_gavage.iterrows():
                file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
                df = pd.read_csv(file_path)            
                df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']
                df_reads_gavage = pd.concat([df_reads_gavage,df])

            #get sample reads dataframe
            #get all samples into dataframe   
            df_samples = df_dose[df_dose['sample_type'] != 'Gavage DNA'].copy()
            #print(df_samples)
            df_sample_reads = pd.DataFrame()
            for idx, row in df_samples.iterrows():
                file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
                df = pd.read_csv(file_path)
                df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']
                df_sample_reads = pd.concat([df_sample_reads,df])

            #pivot all dataframes
            df_sample_pivot = df_sample_reads.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
            df_gavage_pivot = df_reads_gavage.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)

            df_readstable = pd.concat([df_gavage_pivot, df_sample_pivot, df_ntc_pivot], axis = 1).fillna(0)

            df_readstable.columns = df_readstable.columns.str.replace(' ', '-')
            df_readstable.reset_index(inplace = True)

            for strain in li_strains:
                df_readstable.loc[df_readstable['strain'] == strain].drop('strain', axis = 1).to_csv(exportpath / f"{dose}_{strain}_moleculestable.csv", index = False)



    except Exception as e:
        print('Error while creating readtable')
        import traceback
        traceback.print_exc()
    ### P4C3T2
    print(50* '=')
    print(' Creating moleculetables for P4C3T2 10K barcodese')
    print(50*'-')
    print(' ')
    
    basepath = Path(args.basepath)
    exportpath = Path(args.output) / 'P4C3T2/moleculetables'
    exportpath.mkdir(parents=True, exist_ok=True)

    ## Read metadata & ill sheets
    metadir = Path(args.metadir)
    #P4C2T4:
    df_ill = pd.read_excel(metadir / 'ill_sheets' / 'ILL138_GENEWIZ.xls',
                                sheet_name = 'Individual Library Information', header = 1)
    df_samplesheet = pd.read_csv(metadir / 'metadata' / 'P4C3T2_cr_mouse_sample.tsv', sep = '\t')

    df_samplesheet = df_samplesheet.loc[df_samplesheet['sample_type'].str.contains('DNA')].copy()
    df_ill = df_ill.dropna(subset=['Sample Name/Pool Name*'])
    df_ill['sample_id']= df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).split('_')[0])
    df_ill['project'] = df_ill['Library Name*'].apply(lambda x: str(x).split('_')[1])
    #correct library name in ILL138 run
    df_ill['ill'] = df_ill['Library Name*'].apply(lambda x:  'ILL138' if str(x).split('_')[0] == 'ILL139' else str(x).split('_')[0])
    df_ill['file_name'] = df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).replace('_', '-'))
    df_merge = pd.merge(df_ill[['sample_id', 'ill', 'file_name', 'project']], df_samplesheet, how = 'right', on = ['sample_id', 'project'])
    df_merge.head(10)

    df_ill_ntc = df_ill.loc[df_ill['Sample Name/Pool Name*'].str.contains("WATER|NTC")]
    df_merge = df_merge.dropna(subset= ['file_name'])
    print(df_merge)
    print(df_ill_ntc)

    
    li_strains = ['ST1', 'ST2']
    df_ntc = pd.DataFrame()
    for idx, row in df_ill_ntc.iterrows():
        file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
        df = pd.read_csv(file_path)
        df['sample'] = row['Sample Name/Pool Name*']+'_'+row['sample_id']
        df_ntc = pd.concat([df_ntc, df])
    
    df_ntc_pivot = df_ntc.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
    print('FOUND GAVAGE DOSES:')
    print(df_merge['gavage_dose_200ul'].unique())
    for dose in df_merge['gavage_dose_200ul'].unique():
        print(dose)
        df_dose = df_merge.loc[df_merge['gavage_dose_200ul']==dose]
        #get gavage reads dataframe
        #drop files without sequencing data
        df_merge.dropna(subset=['file_name'], inplace = True)
        df_gavage = df_dose[df_dose['sample_type'] == 'Gavage DNA'].copy()
        print(df_gavage)
        df_reads_gavage = pd.DataFrame()
        for idx, row in df_gavage.iterrows():
            file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
            df = pd.read_csv(file_path)            
            df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']
            df_reads_gavage = pd.concat([df_reads_gavage,df])

        #get sample reads dataframe
        #get all samples into dataframe   
        df_samples = df_dose[df_dose['sample_type'] != 'Gavage DNA'].copy()
        #print(df_samples)
        df_sample_reads = pd.DataFrame()
        for idx, row in df_samples.iterrows():
            file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
            df = pd.read_csv(file_path)
            df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']
            df_sample_reads = pd.concat([df_sample_reads,df])

        #pivot all dataframes
        df_sample_pivot = df_sample_reads.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
        df_gavage_pivot = df_reads_gavage.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)

        df_readstable = pd.concat([df_gavage_pivot, df_sample_pivot, df_ntc_pivot], axis = 1).fillna(0)

        df_readstable.columns = df_readstable.columns.str.replace(' ', '-')
        df_readstable.reset_index(inplace = True)

        for strain in li_strains:
            df_readstable.loc[df_readstable['strain'] == strain].drop('strain', axis = 1).to_csv(exportpath / f"{dose}_{strain}_moleculestable.csv", index = False)



    
    ### P4C1T8
    print('Creating Moleculetables for P4C1T8')
    basepath = Path(args.basepath)
    exportpath = Path(args.output) / 'P4C1T8/moleculetables'
    exportpath.mkdir(parents=True, exist_ok=True)

    ## Read metadata & ill sheets
    metadir = Path(args.metadir)
    df_ill = pd.read_excel(metadir / 'ill_sheets' / 'ILL138_GENEWIZ.xls',
                            sheet_name = 'Individual Library Information', header = 1)
    df_samplesheet = pd.read_csv(metadir / 'metadata' / 'P4C1T8-sample_sheets - P4C1T8.tsv', sep = '\t')

    df_samplesheet = df_samplesheet.loc[df_samplesheet['sample_type'].str.contains('DNA')].copy()

    df_ill = df_ill.dropna(subset=['Sample Name/Pool Name*'])
    df_ill['sample_id']= df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).split('_')[0])
    df_ill['project'] = df_ill['Library Name*'].apply(lambda x: str(x).split('_')[1])
    #correct library name in ILL138 run
    df_ill['ill'] = df_ill['Library Name*'].apply(lambda x:  'ILL138' if str(x).split('_')[0] == 'ILL139' else str(x).split('_')[0])
    df_ill['file_name'] = df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).replace('_', '-'))
    df_merge = pd.merge(df_ill[['sample_id', 'ill', 'file_name', 'project']], df_samplesheet, how = 'right', on = ['sample_id', 'project'])

    print(df_merge.head(3))
    print(df_empty.head(3))
    # merge in or map 
    df_merge = pd.merge(df_merge, df_empty, how = 'left', on = ['ill', 'file_name'])
    print(df_merge)
    df_ill_ntc = df_ill.loc[df_ill['Sample Name/Pool Name*'].str.contains("WATER|NTC")]
    
    ## c1 as gavage for cages 1 to 4, cage 5 has own gavage
    strains_c1toc4_mAgavage = ['ST1', 'ST2']
    strains_c1toc4_mBgavage = [ 'ST5', 'ST6']
    strains_c5 = ['ST1', 'ST2', 'ST4']
    #get NTC data
    li_ntc = []
    for idx, row in df_ill_ntc.iterrows():
        file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
        df = pd.read_csv(file_path)
        #manually added empty to negative control
        df['sample'] = row['Sample Name/Pool Name*']+'_'+row['sample_id']+'_empty'
        li_ntc.append(df)
    df_ntc = pd.concat(li_ntc)
    df_ntc_pivot = df_ntc.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
    print(df_ntc_pivot)

    #creates a readtable per mouse gavage, only taking cage1 as gavage
    df_mA_gavage = df_merge[((df_merge['mouse_id'] == 'mA') & (df_merge['sample_type'] == 'gavage DNA')& (df_merge['cage_id'] == 'c1'))]
    df_mB_gavage = df_merge[((df_merge['mouse_id'] == 'mB') & (df_merge['sample_type'] == 'gavage DNA')& (df_merge['cage_id'] == 'c1'))]

    for idx, row in df_mA_gavage.iterrows():
        file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
        df = pd.read_csv(file_path)
        if row['is_empty'] == True:
            df['sample'] =row['project']+'_'+row['file_name']+'_'+row['sample_type']+'_empty'
        else:
            df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']       
        df_reads_mA = df

    for idx, row in df_mB_gavage.iterrows():
        file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
        df = pd.read_csv(file_path) 
        if row['is_empty'] == True:
            df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']+'_empty'
        else:
            df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']       
        df_reads_mB = df

    df_gavage_mA = df_reads_mA.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
    df_gavage_mB = df_reads_mB.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)

    #get cage data
    for cage in df_merge['cage_id'].unique():
        print(cage)
        if cage != 'c5':
            df_cage = df_merge[df_merge['cage_id'] == cage].copy()
            #get all samples into dataframe   
            df_samples = df_cage[df_cage['sample_type'] != 'gavage DNA'].copy()
            #print(df_samples)
            li_reads = []
            for idx, row in df_samples.iterrows():
                file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
                df = pd.read_csv(file_path)
                if row['is_empty'] == True:
                    df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']+'_empty'
                else:
                    df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']    
                li_reads.append(df)            
            df_reads = pd.concat(li_reads)
            
            #pivot all dataframes
            df_reads_pivot = df_reads.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
            
            
            df_reads_mA = pd.concat([df_gavage_mA, df_reads_pivot, df_ntc_pivot], axis = 1).fillna(0)
            df_reads_mB = pd.concat([df_gavage_mB, df_reads_pivot, df_ntc_pivot], axis = 1).fillna(0)
            #replace spaces with '-'
            df_reads_mA.columns = df_reads_mA.columns.str.replace(' ', '-')
            df_reads_mB.columns = df_reads_mB.columns.str.replace(' ', '-')
            df_reads_mA.reset_index(inplace = True)
            df_reads_mB.reset_index(inplace = True)
            
            #export csv files
            for strain in strains_c1toc4_mAgavage:
                df_reads_mA = df_reads_mA.astype({col: 'int' for col in df_reads_mA.select_dtypes(include=['float']).columns})
                output_file = exportpath / f'{cage}_{strain}_moleculestable.csv'
                df_reads_mA.loc[df_reads_mA['strain'] == strain].drop('strain', axis = 1).to_csv(output_file, index = False)
            for strain in strains_c1toc4_mBgavage:
                df_reads_mB = df_reads_mB.astype({col: 'int' for col in df_reads_mB.select_dtypes(include=['float']).columns})
                output_file = exportpath / f'{cage}_{strain}_moleculestable.csv'
                df_reads_mB.loc[df_reads_mB['strain'] == strain].drop('strain', axis = 1).to_csv(output_file, index = False)               
        else:
            df_cage = df_merge[df_merge['cage_id'] == cage].copy()
            
            #creates a readtable per mouse gavage
            df_mA_gavage = df_cage[df_cage['sample_type'] == 'gavage DNA']
            for idx, row in df_mA_gavage.iterrows():
                file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
                df = pd.read_csv(file_path)
                if row['is_empty'] == True:
                    df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']+'_empty'
                else:
                    df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']                
                df_reads_mA = df
            
            #get all samples into dataframe   
            df_samples = df_cage[df_cage['sample_type'] != 'gavage DNA'].copy()
            #print(df_samples)
            li_reads = []
            for idx, row in df_samples.iterrows():
                file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
                df = pd.read_csv(file_path)
                if row['is_empty'] == True:
                    print(row['is_empty'])
                    df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']+'_empty'
                else:
                    df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']
                li_reads.append(df)
            df_reads = pd.concat(li_reads)

            #pivot all dataframes
            df_reads_pivot = df_reads.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
            df_gavage_mA = df_reads_mA.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)

            df_reads_mA = pd.concat([df_gavage_mA, df_reads_pivot, df_ntc_pivot], axis = 1).fillna(0)
            #replace spaces with '-'
            df_reads_mA.columns = df_reads_mA.columns.str.replace(' ', '-')
            df_reads_mA.reset_index(inplace = True)
            
            #export csv files
            for strain in strains_c5:
                df_reads_mA = df_reads_mA.astype({col: 'int' for col in df_reads_mA.select_dtypes(include=['float']).columns})
                output_file = exportpath / f'{cage}_{strain}_moleculestable.csv'
                df_reads_mA.loc[df_reads_mA['strain'] == strain].drop('strain', axis = 1).to_csv(output_file, index = False)

    ### P4C2T4T5
    print('Creating Moleculetables for P4C2T4T5')
    exportpath = Path(args.output) / 'P4C2T4T5/moleculetables'
    exportpath.mkdir(parents=True, exist_ok=True)

    ## P4C2T4:
    df_ill = pd.read_excel(metadir / 'ill_sheets' / 'ILL130_P4C2T4.xlsx',
                                sheet_name = 'Individual Library Information', header = 1)
    df_samplesheet = pd.read_csv(metadir / 'metadata' / 'P4C2T4_cr_mouse_sample - P4C2T4.tsv', sep = '\t')

    df_samplesheet = df_samplesheet.loc[df_samplesheet['sample_type'].str.contains('DNA')].copy()
    df_ill = df_ill.dropna(subset=['Sample Name/Pool Name*'])
    df_ill['sample_id']= df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).split('_')[0])
    df_ill['project'] = df_ill['Library Name*'].apply(lambda x: str(x).split('_')[1])
    df_ill['ill'] = df_ill['Library Name*'].apply(lambda x:  str(x).split('_')[0])
    df_ill['file_name'] = df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).replace('_', '-'))
    df_merge = pd.merge(df_ill[['sample_id', 'ill', 'file_name', 'project']], df_samplesheet, how = 'right', on = ['sample_id', 'project'])
    df_merge.head(10)
    df_P4C2T4 = df_merge.copy()
    df_ill_ntc = df_ill.loc[df_ill['Sample Name/Pool Name*'].str.contains("WATER|NTC")]
    ## P4C2T5
    df_ill = pd.read_excel(metadir / 'ill_sheets' / 'ILL135_sample_sheet.xlsx',
                                sheet_name = 'Individual Library Information', header = 1)
    df_samplesheet = pd.read_csv(metadir / 'metadata' / 'P4C2T5_cr_mouse_sample - P4C2T5.tsv', sep = '\t')

    df_samplesheet = df_samplesheet.loc[df_samplesheet['sample_type'].str.contains('DNA')].copy()
    df_ill = df_ill.dropna(subset=['Sample Name/Pool Name*'])
    df_ill['sample_id']= df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).split('_')[0])
    df_ill['project'] = df_ill['Library Name*'].apply(lambda x: str(x).split('_')[1])
    df_ill['ill'] = df_ill['Library Name*'].apply(lambda x:  str(x).split('_')[0])
    df_ill['file_name'] = df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).replace('_', '-'))
    df_merge = pd.merge(df_ill[['sample_id', 'ill', 'file_name', 'project']], df_samplesheet, how = 'right', on = ['sample_id', 'project'])
    df_P4C2T5 = df_merge.copy()
    df_ill_ntc = pd.concat([df_ill_ntc, df_ill.loc[df_ill['Sample Name/Pool Name*'].str.contains("WATER|NTC")]], axis= 0)

    #merge
    df_merge = pd.concat([df_P4C2T4, df_P4C2T5], axis = 0)

    li_strains = ['ST1', 'ST2', 'ST4', 'ST5', 'ST6']

    df_ntc = pd.DataFrame()
    for idx, row in df_ill_ntc.iterrows():
        file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
        df = pd.read_csv(file_path)
        df['sample'] = row['Sample Name/Pool Name*']+'_'+row['sample_id']
        df_ntc = pd.concat([df_ntc, df])
    df_ntc_pivot = df_ntc.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)


    #get gavage reads dataframe
    #drop files without sequencing data
    df_merge.dropna(subset=['file_name'], inplace = True)
    df_gavage = df_merge[df_merge['sample_type'] == 'Gavage DNA'].copy()
    df_reads_gavage = pd.DataFrame()
    for idx, row in df_gavage.iterrows():
        file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
        df = pd.read_csv(file_path)
        df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']
        df_reads_gavage = pd.concat([df_reads_gavage,df])

    #get sample reads dataframe
    df_samples = df_merge[df_merge['sample_type'] != 'Gavage DNA'].copy()
    #print(df_samples)
    df_sample_reads = pd.DataFrame()
    for idx, row in df_samples.iterrows():
        file_path = basepath / row['ill'] / 'readtables' / row['file_name'] / f"{row['file_name']}_total_moleculestable.csv"
        df = pd.read_csv(file_path)
        df['sample'] = row['project']+'_'+row['file_name']+'_'+row['sample_type']
        df_sample_reads = pd.concat([df_sample_reads,df])

    #pivot all dataframes
    df_sample_pivot = df_sample_reads.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
    df_gavage_pivot = df_reads_gavage.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)

    df_readstable = pd.concat([df_gavage_pivot, df_sample_pivot, df_ntc_pivot], axis = 1).fillna(0)

    df_readstable.columns = df_readstable.columns.str.replace(' ', '-')
    df_readstable.reset_index(inplace = True)

    for strain in li_strains:
        df_readstable = df_readstable.astype({col: 'int' for col in df_readstable.select_dtypes(include=['float']).columns})
        output_file = exportpath / f'{strain}_moleculestable.csv'
        df_readstable.loc[df_readstable['strain'] == strain].drop('strain', axis = 1).to_csv(output_file, index = False)


        ### P4C4T1
    print('Creating Moleculetabel for P4C4T1 In Vitro')
    basepath = Path(args.basepath)
    exportpath = Path(args.output) / 'P4C4T1/moleculetables'
    exportpath.mkdir(parents=True, exist_ok=True)
    metadir = Path(args.metadir)
    df_ill = pd.read_excel(metadir / 'ill_sheets' / 'ILL132_P4C4T1.xlsx',
                            sheet_name = 'Individual Library Information', header = 1)

    df_samplesheet = pd.read_csv(metadir / 'metadata' / 'P4C4T1_metadata.tsv', sep = '\t')
    df_samplesheet['project'] = 'P4C4T1'
    df_ill = df_ill.dropna(subset=['Sample Name/Pool Name*'])
    df_ill['sample_id']= df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).split('_')[0])
    df_ill['project'] = df_ill['Library Name*'].apply(lambda x: str(x).split('_')[1])
    #correct library name in ILL138 run
    df_ill['ill'] = df_ill['Library Name*'].apply(lambda x:  'ILL138' if str(x).split('_')[0] == 'ILL139' else str(x).split('_')[0])
    df_ill['file_name'] = df_ill['Sample Name/Pool Name*'].apply(lambda x: str(x).replace('_', '-'))
    df_merge = pd.merge(df_ill[['sample_id', 'ill', 'file_name', 'project']], df_samplesheet, how = 'right', on = ['sample_id', 'project'])
    df_merge.head(10)

    df_ill_ntc = df_ill.loc[df_ill['Sample Name/Pool Name*'].str.contains("WATER|NTC")]
    df_merge = df_merge.dropna(subset= ['file_name'])
    print(df_merge)
    print(df_ill_ntc)

    # generate readtables for getFP
    li_strains = ['ST1', 'ST2', 'ST4', 'ST5', 'ST6']
    #get gavage reads dataframe
    #drop files without sequencing data
    df_merge.dropna(subset=['file_name'], inplace = True)

    #print(df_samples)
    df_sample_reads = pd.DataFrame()
    for idx, row in df_merge.iterrows():
        df = pd.read_csv(basepath /  row['ill'] / 'readtables' /  row['file_name'] / f"{row['file_name']}_total_moleculestable.csv")
        df['sample'] = row['project']+'_'+row['file_name']
        df_sample_reads = pd.concat([df_sample_reads,df])

    #pivot all dataframes
    df_sample_pivot = df_sample_reads.pivot(index = ['strain', 'umi_seq'], columns = 'sample', values = 'molecules').fillna(0)
    df_readstable = df_sample_pivot.copy()
    df_readstable.columns = df_readstable.columns.str.replace(' ', '-')
    df_readstable.reset_index(inplace = True)
    for strain in li_strains:
        df_readstable.loc[df_readstable['strain'] == strain].drop('strain', axis = 1).to_csv(exportpath / f"{strain}_moleculestable.csv", index = False)

    # read all readtables, add gavage columns and drop rows in which all gavage columns are 0
    print('Adding gavage from P4C2T4T5')
    li_strains = ['ST1', 'ST2', 'ST4', 'ST5', 'ST6']
    gavagepath = Path(args.output) / 'P4C2T4T5/moleculetables'
    for strain in li_strains:
        df_gavage = pd.read_csv(gavagepath / f"{strain}_moleculestable.csv")
        #get gavage columns
        li_gavage = [col for col in df_gavage.columns if 'gavage' in col.lower()]
        print(df_gavage[li_gavage])
        df = pd.read_csv(exportpath / f"{strain}_moleculestable.csv")
        df.set_index('umi_seq', inplace = True)
        df_gavage.set_index('umi_seq', inplace = True)
        df = pd.concat([df_gavage[li_gavage] ,df], axis = 1)
        df.reset_index(inplace = True)
        df.fillna(0, inplace = True)
        #drop rows in which all gavage columns are 0
        df[df[li_gavage].sum(axis = 1) > 0].to_csv(exportpath / f"{strain}_moleculestable.csv", index = False)
        
if __name__ == "__main__":
    main()