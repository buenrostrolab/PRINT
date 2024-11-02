#!/bin/bash
#SBATCH -n 16                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-8:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test   # Partition to submit to
#SBATCH --mem=80G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o run_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e run_%j.err  # File to which STDERR will be written, %j inserts jobid

module load samtools/1.10-fasrc01
module load bedtools2/2.26.0-fasrc01

cores=16
cd /n/holylfs05/LABS/buenrostro_lab/Users/yanhu/cell_dynamics/multiScaleFootprinting/data/humanGenomicDNA/

if [ ! -d fragments ]; then
  mkdir fragments
fi

cd bams

bam_files=*.st.rmdup.flt.bam
for bam in $bam_files
do

    fname=`echo $bam | sed 's/_S[0-9]_001.st.rmdup.flt.bam//g'`
    echo "Processing $fname"
    
    if [ ! -f $fname.n_sorted.rmdup.flt.bam ]; then
      echo "Name sorting bam file"
      samtools sort -n -@ $cores -o $fname.n_sorted.rmdup.flt.bam $fname*.st.rmdup.flt.bam
    fi
  
    if [ ! -f ../fragments/$fname.frags.gz ]; then
      echo "Converting bam file to fragments file"
      bedtools bamtobed -i $fname.n_sorted.rmdup.flt.bam -bedpe | \
      sed 's/_/\t/g' | \
      awk -v OFS="\t" -v ID=$fname '{if($9=="+"){print $1,$2+4,$6-5,ID}}' |\
      sort --parallel=$cores -S 40G  -k4,4 -k1,1 -k2,2n -k3,3n | \
      uniq -c | \
      awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | \
      gzip > $fname.frags.gz
      mv $fname.frags.gz ../fragments
    fi
  
done