#!/bin/bash
#SBATCH -t 4-0:0:0
#SBATCH -c 1
module load snakemake

#snakemake dedupl_UMI --snakefile Snakefile \
#  --configfile config/general.yaml \
#  --cluster-config config/config.json \
#  --cluster "sbatch -p {cluster.partition} -c {threads} -t {cluster.time} --mem-per-cpu={cluster.mem_per_cpu}" \
#  --jobs 50 \
#  --rerun-incomplete \
#  --latency-wait 60

#snakemake merge_hs_metrics_UMI_start_stop --snakefile Snakefile \
#  --configfile config/general.yaml \
#  --cluster-config config/config.json \
#  --cluster "sbatch -p {cluster.partition} -c {threads} -t {cluster.time} --mem-per-cpu={cluster.mem_per_cpu}" \
#  --jobs 50 \
#  --rerun-incomplete \
#  --latency-wait 60

#snakemake picard_QC_UMI_start_stop --snakefile Snakefile \
#  --configfile config/general.yaml \
#  --cluster-config config/config.json \
#  --cluster "sbatch -p {cluster.partition} -c {threads} -t {cluster.time} --mem-per-cpu={cluster.mem_per_cpu}" \
#  --jobs 50 \
#  --rerun-incomplete \
#  --latency-wait 60

#snakemake get_consensus_mapped_filtered_bam --snakefile Snakefile \
#  --configfile config/general.yaml \
#  --cluster-config config/config.json \
#  --cluster "sbatch -p himem -c {threads} -t {cluster.time} --mem-per-cpu={cluster.mem_per_cpu}" \
#  --jobs 50 \
#  --rerun-incomplete \
#  --latency-wait 60

#snakemake get_Mutect2_tumourOnly --snakefile Snakefile \
#  --configfile config/general.yaml \
#  --cluster-config config/config.json \
#  --cluster "sbatch -p all -c {threads} -t {cluster.time} --mem-per-cpu=30G" \
#  --jobs 50 \
#  --rerun-incomplete \
#  --latency-wait 60

#snakemake get_Mutect2_filtered --snakefile Snakefile \
#  --configfile config/general.yaml \
#  --cluster-config config/config.json \
#  --cluster "sbatch -p all -c {threads} -t {cluster.time} --mem-per-cpu=30G" \
#  --jobs 50 \
#  --rerun-incomplete \
#  --latency-wait 60

#snakemake get_Funcotator --snakefile Snakefile \
#  --configfile config/general.yaml \
#  --cluster-config config/config.json \
#  --cluster "sbatch -p all -c {threads} -t {cluster.time} --mem-per-cpu=30G" \
#  --jobs 50 \
#  --rerun-incomplete \
#  --latency-wait 60

#snakemake get_CNVs --snakefile Snakefile \
#  --configfile config/general.yaml \
#  --cluster-config config/config.json \
#  --cluster "sbatch -p all -c {threads} -t {cluster.time} --mem=30G" \
#  --jobs 50 \
#  --rerun-incomplete \
#  --latency-wait 60

#snakemake cnvkit_heatmap --snakefile Snakefile \
#  --configfile config/general.yaml \
#  --cluster-config config/config.json \
#  --cluster "sbatch -p all -c {threads} -t {cluster.time} --mem=30G" \
#  --jobs 50 \
#  --rerun-incomplete \
#  --latency-wait 60

snakemake cnvkit --snakefile Snakefile \
  --configfile config/general.yaml \
  --cluster-config config/config.json \
  --cluster "sbatch -p all -c {threads} -t {cluster.time} --mem=30G" \
  --jobs 50 \
  --rerun-incomplete \
  --latency-wait 60
