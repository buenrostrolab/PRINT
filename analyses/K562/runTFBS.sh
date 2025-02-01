#!/bin/bash
#SBATCH -n 16                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 2-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p bigmem   # Partition to submit to
#SBATCH --mem=300G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o run_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e run_%j.err  # File to which STDERR will be written, %j inserts jobid

module load R/4.1.0-fasrc01
singularity exec --cleanenv --env R_LIBS_USER=/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/apps/R_4.1.0 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_13.sif Rscript main.R 1
