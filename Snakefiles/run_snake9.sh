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

#remove old logs
#rm -rf /sc/arion/work/ruprec01/results/P4C2-run1/FP
#rm -rf /sc/arion/work/ruprec01/results/P4C2-run1/FP_norm
#rm -rf /sc/arion/work/ruprec01/results/P4C2-run1/readtables
#rm -rf /sc/arion/work/ruprec01/results/P4C2-run1/heatmaps
#rm -rf /sc/arion/work/ruprec01/results/P4C2-run1/distance
#rm -rf /sc/arion/work/ruprec01/results/P4C2-run1-merge

# remove clustger log files
rm -rf /sc/arion/work/ruprec01/log/cluster/*
mkdir -p /sc/arion/work/ruprec01/log/cluster


#conda activate snakemake

# Run with Snakemake 9.1.1
# rerung specific rules -R get_single_GD get_GD_aggregate
snakemake \
  -s ${SNAKEFILE} \
  --jobs 200 \
  --cores 1 \
  --rerun-incomplete \
  --latency-wait 240 \
  --executor lsf \
  --default-resources \
    mem_mb=64000 \
    disk_mb=30000 \
    lsf_project="acc_faithj02a" \
    lsf_queue="premium" \
    walltime=720 \
    "lsf_extra='-o /sc/arion/work/ruprec01/log/cluster/%J.out -e /sc/arion/work/ruprec01/log/cluster/%J.err -L /bin/bash'"