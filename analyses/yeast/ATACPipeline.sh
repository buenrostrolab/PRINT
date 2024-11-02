#!/bin/bash
#SBATCH -n 12                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-08:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p bigmem   # Partition to submit to
#SBATCH --mem=300G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o run_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e run_%j.err  # File to which STDERR will be written, %j inserts jobid

# Load required modules
module load R/4.1.0-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.1.0:$R_LIBS_USER
module load macs2/2.1.1.20160309-fasrc02
module load bowtie2/2.3.2-fasrc02
module load samtools/1.10-fasrc01
module load bedtools2/2.26.0-fasrc01
cores=16
cd /n/holylfs05/LABS/buenrostro_lab/Users/yanhu/cell_dynamics/multiScaleFootprinting/data/yeast/ATAC/bam

for ID in {1822137..1822147}
do

  fname=yeast${ID}_S1
  
  if [ ! -f $fname.n_sorted.rmdup.flt.bam ]; then
    echo "Name sorting bam file $fname.sorted.bam"
    samtools sort -n -@ $cores -o $fname.n_sorted.rmdup.flt.bam $fname.st.rmdup.flt.bam 
  fi
    
  if [ ! -f yeast${ID}.frags.gz ]; then
    echo "Converting bam file to fragments file"
    bedtools bamtobed -i $fname.n_sorted.rmdup.flt.bam -bedpe | \
    sed 's/_/\t/g' | \
    awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6-5}}' |\
    sort --parallel=$cores -S 40G -k1,1 -k2,2n -k3,3n | \
    uniq -c | \
    awk -v OFS="\t" '{print $2, $3, $4, $1}' | \
    gzip > yeast${ID}.frags.gz
  fi
    
done

# Merge fragments files
if [ ! -f merged.frags.gz ]; then
  zcat *.frags.gz | gzip > merged.frags.gz
fi

# Peak Calling
if [ ! -f peaks.bed ]; then

  macs2 callpeak --nomodel -t merged.frags.gz \
        --outdir ./peakCalling -n yeast -f BEDPE \
        --nolambda --keep-dup all --call-summits
        
fi

# Re-organize folders
mv merged.frags.gz ../../
mv peakCalling ../
mkdir ../frags
mv *frags* ../frags