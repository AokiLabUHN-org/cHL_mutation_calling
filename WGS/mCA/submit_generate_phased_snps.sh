#!/bin/bash
#SBATCH -t 3-0:0:0
#SBATCH -p veryhimem
#SBATCH --mem 160G
#SBATCH --cpus-per-task=8

module load parallel
module load samtools
module load bwa

cd ~/mCAs_WGS/pipeline
bash generate_phased_snps.sh \
  -m cHL_WGS.test.manifest.txt \
  -o /cluster/projects/aokigroup/Stueckmann/out/phased_snps_test \
  -n cHL_WGS_test \
  --remap -j 8 -J 1 \
  -c chr21 -v
