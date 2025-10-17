#!/bin/bash
# Setup

# Default Snakefile
SNAKEFILE="Snakefiles/ampliconez_Snakefile"

# Parse command-line arguments
while getopts "s:" opt; do
  case $opt in
    s)
      SNAKEFILE="$OPTARG"
      ;;
    *)
      echo "Usage: $0 [-s snakefile_path]"
      exit 1
      ;;
  esac
done

# remove clustger log files
rm -rf /sc/arion/work/ruprec01/log/cluster/*
mkdir -p /sc/arion/work/ruprec01/log/cluster


#Define inout and output directories
INPUT_DIRS=('/sc/arion/projects/faithj02a/data/seq/2025/ILL138'
            '/sc/arion/projects/faithj02a/data/seq/2025/ILL135'
            '/sc/arion/projects/faithj02a/data/seq/2025/ILL134'
            '/sc/arion/projects/faithj02a/data/seq/2025/ILL132/30-1179742904/00_fastq'
             '/sc/arion/projects/faithj02a/data/seq/2025/ILL130/00_fastq'
            )
OUTPUT_DIRS=('/sc/arion/work/ruprec02/population-bottlenecks/ILL138'
             '/sc/arion/work/ruprec02/population-bottlenecks/ILL135'
             '/sc/arion/work/ruprec02/population-bottlenecks/ILL134'
             '/sc/arion/work/ruprec02/population-bottlenecks/ILL132'
             '/sc/arion/work/ruprec02/population-bottlenecks/ILL130')
SCRATCH_DIR=('/sc/arion/scratch/ruprec02/population-bottlenecks/ILL138'
             '/sc/arion/scratch/ruprec02/population-bottlenecks/ILL135'
             '/sc/arion/scratch/ruprec02/population-bottlenecks/ILL134'
             '/sc/arion/scratch/ruprec02/population-bottlenecks/ILL132'
             '/sc/arion/scratch/ruprec02/population-bottlenecks/ILL130'
             )

# Run with Snakemake 9.1.1
# rerung specific rules -R get_single_GD get_GD_aggregate
for i in "${!INPUT_DIRS[@]}"; do
  snakemake \
    -s ${SNAKEFILE} \
    --config input_dir="${INPUT_DIRS[$i]}" output_dir="${OUTPUT_DIRS[$i]}" scratch_dir="${SCRATCH_DIR[$i]}" \
    --jobs 500 \
    --cores 1 \
    --rerun-incomplete \
    --latency-wait 240 \
    --executor lsf \
    --default-resources \
      mem_mb=40000 \
      disk_mb=30000 \
      lsf_project="acc_faithj02a" \
      lsf_queue="express" \
      walltime=480 \
      "lsf_extra='-o /sc/arion/work/ruprec01/log/cluster/%J.out -e /sc/arion/work/ruprec01/log/cluster/%J.err -L /bin/bash'"
done