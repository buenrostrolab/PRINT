#!/bin/bash
#SBATCH -n 12                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 1-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p bigmem   # Partition to submit to
#SBATCH --mem=200G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o run_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e run_%j.err  # File to which STDERR will be written, %j inserts jobid

cd /n/home09/yanhu/cell_dynamics/multiScaleFootprinting/data/yeast/nucleoATAC

nucleoatac run \
  --bed ../peaks.bed \
  --bam yeast.merged.st.rmdup.flt.bam \
  --fasta /n/home09/yanhu/cell_dynamics/multiScaleFootprinting/data/shared/refGenomes/sacCer3.fa \
  --out yeastCRE