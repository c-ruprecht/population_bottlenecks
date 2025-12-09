import os
from snakemake.io import expand, directory
import re
import glob
import pandas as pd

# PATHS
OUTPUT_DIR = '/sc/arion/work/ruprec02/population-bottlenecks/'
SCRIPTS_DIR = '/hpc/users/ruprec02/git/population_bottlenecks/scripts'
METADIR = '/hpc/users/ruprec02/git/population_bottlenecks/input_files'

# Create main output directory structure
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load FP configuration from TSV
FP_CONFIG_FILE = os.path.join(METADIR, 'nb_config.tsv')
fp_config_df = pd.read_csv(FP_CONFIG_FILE, sep='\t')

# Create lookup dictionary: (experiment, strain, anal_folder) -> params
fp_config = {}
for _, row in fp_config_df.iterrows():
    key = (row['experiment'], row['strain'], row['anal_folder'])
    fp_config[key] = {
        'umi_threshold': int(row['umi_threshold'])
    }

print(f"Loaded {len(fp_config)} FP configurations from {FP_CONFIG_FILE}")


rule all:
    input:
        # Nb analysis outputs based on config file
        expand(OUTPUT_DIR + "FP/{experiment}/{anal_folder}/{strain}_Nb_estimates.csv",
               zip,
               experiment=fp_config_df['experiment'].tolist(),
               strain=fp_config_df['strain'].tolist(),
               anal_folder=fp_config_df['anal_folder'].tolist())

checkpoint create_experiments:
    output:
        directory(OUTPUT_DIR + "moleculetables")
    params:
        basepath = OUTPUT_DIR,
        metadir = METADIR,
        scripts_dir = SCRIPTS_DIR,
        empty_samples= METADIR+'/non_empty_sample.csv'
    log:
        OUTPUT_DIR + "logs/create_moleculetables.log"
    shell:
        """
        conda run -n umi_tools_env python {params.scripts_dir}/create_experiment_moleculetables.py \
            -b {params.basepath} \
            -m {params.metadir} \
            -o {output} \
            -e {params.empty_samples} &> {log}
        """

# defines input molecules tables from the config file for moleculetable clustering
def get_moleculetable_input(wildcards):
    """Input function that waits for checkpoint completion"""
    checkpoints.create_experiments.get()
    return OUTPUT_DIR + f"moleculetables/{wildcards.experiment}/moleculetables/{wildcards.strain}_moleculestable.csv"

rule cluster_moleculetables:
    input:
        readstable = get_moleculetable_input
    output:
        OUTPUT_DIR + "clustered/{experiment}/{anal_folder}/{strain}_clustered.csv",
    params:
        scripts_dir = SCRIPTS_DIR,
        output_prefix = OUTPUT_DIR + "clustered/{experiment}/{anal_folder}/{strain}",
        threshold = lambda wildcards: fp_config[(wildcards.experiment, wildcards.strain, wildcards.anal_folder)]['umi_threshold'],
    log:
        OUTPUT_DIR + "logs/{experiment}/{anal_folder}/{strain}_clustering.log"
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        conda run -n umi_tools_env python {params.scripts_dir}/cluster_moleculetables.py \
            --input {input.readstable} \
            --output {params.output_prefix} \
            --threshold {params.threshold} &> {log}
        """


# Run the Nb calculation
rule get_FP:
    input:
        readstable = OUTPUT_DIR + "clustered/{experiment}/{anal_folder}/{strain}_clustered.csv"
    output:
        OUTPUT_DIR + "FP/{experiment}/{anal_folder}/{strain}_Nb_estimates.csv"
    params:
        output_dir = lambda wildcards: OUTPUT_DIR + f"FP/{wildcards.experiment}/{wildcards.anal_folder}/{wildcards.strain}",
        #output_dir = lambda wildcards: fp_config[(wildcards.experiment, wildcards.strain, wildcards.anal_folder)]['anal'],
        scripts_dir = SCRIPTS_DIR
    log:
        OUTPUT_DIR + "logs/{experiment}/{anal_folder}/{strain}_Nb.log"
    shell:
        """
        # Initialize log file with timestamp
        echo "=== Nb Analysis Started: $(date) ===" > {log}
        echo "Input file: {input.readstable}" >> {log}
        echo "Output directory: {params.output_dir}" >> {log}
        echo "Analysis folder: {wildcards.anal_folder}" >> {log}

        # Create directories
        echo "Creating output directories..." >> {log}
        mkdir -p $(dirname {log}) 2>> {log}
        echo "Directories created successfully" >> {log}

        # Check if input file exists
        echo "Checking input file..." >> {log}
        if [ ! -f "{input.readstable}" ]; then
            echo "ERROR: Input file does not exist: {input.readstable}" >> {log}
            exit 1
        fi
        echo "Input file exists" >> {log}

        # Dynamically find gavage column indices (0-indexed for Python)
        echo "Finding gavage columns..." >> {log}
        WHERE_REFS=$(conda run -n umi_tools_env python -c "
import pandas as pd
import sys
try:
    df = pd.read_csv('{input.readstable}', index_col=0)
    print('Columns in dataframe:', list(df.columns), file=sys.stderr)
    gavage_cols = [i for i, col in enumerate(df.columns) if 'gavage' in str(col).lower()]
    print('Found gavage columns at indices:', gavage_cols, file=sys.stderr)
    print(','.join(map(str, gavage_cols)) if gavage_cols else '0')
except Exception as e:
    print('ERROR finding gavage columns:', str(e), file=sys.stderr)
    print('0')
" 2>> {log})

        echo "Gavage column indices (0-based for Python): $WHERE_REFS" >> {log}

        # Run Nb calculation
        conda run -n umi_tools_env python {params.scripts_dir}/getFP_Nb.py \
            {input.readstable} \
            "$WHERE_REFS" \
            {params.output_dir} 2>&1 | tee -a {log}

        EXIT_CODE=${{PIPESTATUS[0]}}

        # Check script completion
        if [ $EXIT_CODE -eq 0 ]; then
            echo "=== Python script completed successfully: $(date) ===" >> {log}
            echo "Output file created at {output}" >> {log}
        else
            echo "=== Python script FAILED with exit code $EXIT_CODE: $(date) ===" >> {log}
            exit $EXIT_CODE
        fi
        """