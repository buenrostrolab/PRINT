#!/bin/bash
#SBATCH -n 16                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-8:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test   # Partition to submit to
#SBATCH --mem=80G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o run_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e run_%j.err  # File to which STDERR will be written, %j inserts jobid

dataset=HepG2
cd /n/holylfs05/LABS/buenrostro_lab/Users/yanhu/cell_dynamics/multiScaleFootprinting/data/HepG2/

# Run ATACorrect
TOBIAS ATACorrect \
  --bam HepG2.bam \
  --genome ../shared/refGenomes/hg38.fa \
  --peaks peaks.bed \
  --outdir Tobias \
  --cores 16
  
# Run FootprintScores
TOBIAS FootprintScores \
  --signal Tobias/HepG2_corrected.bw \
  --regions peaks.bed \
  --output HepG2_footprints.bw \
  --cores 16
  
mv HepG2_footprints.bw Tobias/

# Run BINDetect
TOBIAS BINDetect \
  --motifs ../shared/cisBP2021_human.jaspar \
  --signals Tobias/HepG2_footprints.bw \
  --genome ../shared/refGenomes/hg38.fa \
  --peaks peaks.bed \
  --outdir Tobias/prediction \
  --cond_names HepG2 --cores 16
