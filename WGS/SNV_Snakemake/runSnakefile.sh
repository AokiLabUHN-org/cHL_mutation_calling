#!/bin/bash
#SBATCH -t 4-0:0:0
#SBATCH -c 1
module load snakemake

#snakemake get_bqsr --snakefile Snakefile \
#  --configfile general.yaml \
#  --cluster-config config.json \
#  --cluster "sbatch -p himem -c 1 -t {cluster.time} --mem=30G" \
#  --jobs 53 \
#  --ignore-incomplete \
#  --latency-wait 60

snakemake get_Funcotator --snakefile Snakefile \
  --configfile general.yaml \
  --cluster-config config.json \
  --cluster "sbatch -p himem -c 1 -t 7-0:0:0 --mem=60G" \
  --jobs 53 \
  --ignore-incomplete \
  --latency-wait 60
