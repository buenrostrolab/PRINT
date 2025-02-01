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
module load macs2/2.1.1.20160309-fasrc02
module load intel/2017c impi/2017.4.239 Bowtie2/2.3.4.1
module load Java/1.8

cores=16
experiment=IFNStim
genome=mm10
ref=/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/utils/bowtie2/genomes/mm10/mm10
picardPATH=/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/utils/picard/picard.jar
cd /n/holylfs05/LABS/buenrostro_lab/Users/yanhu/cell_dynamics/multiScaleFootprinting/data/IFNStim/

if [ ! -d fragments ]; then
  mkdir fragments
fi

cd bams

SRR_files=SRR293909?
for fname in $SRR_files
do

  # Convert SRR files to fastq files
  echo "Converting SRR files to fastq files"
  if [ ! -f $fname.fastq ]; then
    fasterq-dump $fname -e 16
  fi
  
  # Trim adapters
  echo "Trimming adapters"
  if [ ! -f $fname.trimmed.fastq ]; then
    cutadapt -a CTGTCTCTTATACACATCT -o $fname.trimmed.fastq $fname.fastq
  fi
  
  # Align to genome
  echo "Aligning to genome"
  if [ ! -f $fname.sam ]; then
    bowtie2 -X 2000 -p $cores --rg-id $fname -x $ref -U $fname.trimmed.fastq -S $fname.sam
  fi
  
  # Convert sam to bam
  echo "Converting sam file to bam file"
  if [ ! -f $fname.bam ]; then
    samtools view -@ $cores -h -S -b -o $fname.bam $fname.sam
  fi
  
  # Sort bam file
  echo "Sorting bam file"
  if [ ! -f $fname.sorted.bam ]; then
    samtools sort -@ $cores -o $fname.sorted.bam $fname.bam
  fi
  
  # Index bam file
  echo "Indexing the sorted bam file"
  if [ ! -f $fname.sorted.bam.bai ]; then
    samtools index -@ $cores $fname.sorted.bam 
  fi
  
  # Remove dups and mito
  echo "Removing duplicates and unwanted chrs: "
  if [ ! -f $fname.sorted.rmdup.flt.bam ]; then
    chrs=`samtools view -H $fname.sorted.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v Y | awk '{if(length($0)<6)print}'`
    # If paired end, use samtools view -@ $cores -b -q 30 -f 0x2 $fname.sorted.bam -o temp.bam `echo $chrs`
    samtools view -@ $cores -b -q 30 -F 260 $fname.sorted.bam -o temp.bam `echo $chrs`
    java -jar -Djava.io.tmpdir=`pwd`/tmp $picardPATH \
      MarkDuplicates \
      INPUT=temp.bam \
      OUTPUT=$fname.sorted.rmdup.flt.bam \
      METRICS_FILE=$fname.dups.log \
      REMOVE_DUPLICATES=true \
      VALIDATION_STRINGENCY=SILENT
    rm temp.bam
  fi
  
  # Index bam file
  echo "Indexing the sorted bam file"
  if [ ! -f $fname.sorted.rmdup.flt.bam.bai ]; then
    samtools index $fname.sorted.rmdup.flt.bam -@ $cores
  fi
  
  # Name sorting bam file
  echo "Name sorting bam file"
  if [ ! -f $fname.n_sorted.rmdup.flt.bam ]; then
    samtools sort -n -@ $cores -o $fname.n_sorted.rmdup.flt.bam $fname.sorted.rmdup.flt.bam
  fi

  # Convert bam file to fragments file
  echo "Converting bam file to fragments file"
  if [ ! -f ../fragments/$fname.frags.gz ]; then
    bedtools bamtobed -i $fname.n_sorted.rmdup.flt.bam | \
    sed 's/_/\t/g' | \
    awk -v OFS="\t" -v ID=$fname '{if($6=="+"){print $1,$2+4,$2+4,ID}else if($6=="-"){print $1,$3-5,$3-5,ID}}' |\
    sort --parallel=$cores -S 40G  -k4,4 -k1,1 -k2,2n -k3,3n | \
    gzip > $fname.frags.gz
    mv $fname.frags.gz ../fragments
  fi
  
done

# Peak Calling
if [ ! -f ../peaks.bed ]; then

  # Call peaks using MACS2
  mkdir peakCalling
  macs2 callpeak \
    -t *.sorted.rmdup.flt.bam \
    -f BAM \
    -n $experiment \
    --outdir peakCalling \
    --keep-dup all \
    --nolambda --nomodel \
    --call-summits 
    
  mv peakCalling ../
  cd ../
  
  # Filter and resize peaks
  singularity exec --cleanenv --env R_LIBS_USER=/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/apps/R_4.1.0 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_13.sif Rscript ../../code/filterPeaks.R $genome
  
fi

# Merge fragments files
cd fragments
if [ ! -f ../all.frags.tsv.gz ]; then
  zcat *frags.gz | gzip > all.frags.tsv.gz
  mv all.frags.tsv.gz ../
fi