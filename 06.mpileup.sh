#!/usr/bin/env bash 
# Mpileup SNP Caller Script Creator
# Chad Smith
#
# Outputs compressed VCF file of SNPs from BAM alignments 
# -Searches for bam files in $BAM_DIR and divides files up for the cluster. How files are divided up
# is determined by the user. Samtools uses population allele frequencies to call SNPs, so a higher
# number of samples will increase the effectiveness of SNP calling at the cost of longer run times
# on the cluster
# -Creates script to make *.bai index files
# -Runs samtools mpileup to call SNPs
#
# BASE		working directory
# BAM_DIR	dir of bam files to be analyzed
# OUT_DIR	output directory 
# GENOME	Path to genome
# SAMTOOLS	path to samtools
# BCFTOOLS	path to bcftools
#
# launcher_creator.py must be in $PATH. This creates a script for TACC to launch the jobs. 
#
# Usage:
# ./06.mpileup.sh

BASE=/scratch/02716/cs762/atex_popgen
OUT_DIR=${BASE}/mpileup/q20/
BAM_DIR=${BASE}/bwa/q20/
GENOME=/work/02716/cs762/genomes/atex_v0.1.fa
SAMTOOLS=/home1/02716/cs762/.local/bin/samtools
BCFTOOLS=/home1/02716/cs762/.local/bin/bcftools

if [ ! -d $OUT_DIR ]; then mkdir -p $OUT_DIR;fi

# Collate names of fq files and divides them for the cluster.
find ${BAM_DIR} -maxdepth 1 -iname "*.bam" > /tmp/06.bam_all.txt
sed -n '1,87p' /tmp/06.bam_all.txt > 06.bam_1.txt
sed -n '88,174p' /tmp/06.bam_all.txt > 06.bam_2.txt
sed -n '175,261p' /tmp/06.bam_all.txt > 06.bam_3.txt
sed -n '262,348p' /tmp/06.bam_all.txt > 06.bam_4.txt
sed -n '349,432p' /tmp/06.bam_all.txt > 06.bam_5.txt

# Outputs a script to create BAM bai index files.
if [ -f bai.sh ]; then rm bai.sh; fi

while read line; do 
 if [ ! -f ${line}.bai ]; then 
  echo "$SAMTOOLS index $line" >> bai.sh
 fi
done < /tmp/06.bam_all.txt

# TACC script for indexing BAM files
if [ -f bai.sh ];then launcher_creator.py -n bai -t 02:00:00 -q normal -w 24 -j bai.sh -l bai.slurm;fi

if [ -f 06.mpileup_ex.sh ];then rm 06.mpileup_ex.sh;fi
nfq=`wc -l /tmp/06.bam_all.txt | cut -f1 -d " "`
echo "$nfq samples to process" 

# $SAMTOOLS  -C Coefficient for downgrading mapping quality - 50 recommended for bwa
#               -d maximum depth at any position
#               -BI Reduce false SNPs from misalignments, no Indels
#               -g Genotype and output to BCF
# $BCFTOOLS -vc -O z output variant sites only, use consensus caller, output=compressed vcf
echo "$SAMTOOLS mpileup -C 50 -d 500 -t AD,INFO/AD,DP -BI -g -f $GENOME -b 06.bam_1.txt |\
$BCFTOOLS call -vc -O z -f GQ -o $OUT_DIR/atex_raw1.vcf.gz" >> 06.mpileup_ex.sh
echo "$SAMTOOLS mpileup -C 50 -d 500 -t AD,INFO/AD,DP -BI -g -f $GENOME -b 06.bam_2.txt |\
$BCFTOOLS call -vc -O z -f GQ -o $OUT_DIR/atex_raw2.vcf.gz" >> 06.mpileup_ex.sh
echo "$SAMTOOLS mpileup -C 50 -d 500 -t AD,INFO/AD,DP -BI -g -f $GENOME -b 06.bam_3.txt |\
$BCFTOOLS call -vc -O z -f GQ -o $OUT_DIR/atex_raw3.vcf.gz" >> 06.mpileup_ex.sh
echo "$SAMTOOLS mpileup -C 50 -d 500 -t AD,INFO/AD,DP -BI -g -f $GENOME -b 06.bam_4.txt |\
$BCFTOOLS call -vc -O z -f GQ -o $OUT_DIR/atex_raw4.vcf.gz" >> 06.mpileup_ex.sh
echo "$SAMTOOLS mpileup -C 50 -d 500 -t AD,INFO/AD,DP -BI -g -f $GENOME -b 06.bam_5.txt |\
$BCFTOOLS call -vc -O z -f GQ -o $OUT_DIR/atex_raw5.vcf.gz" >> 06.mpileup_ex.sh

# TACC script for calling SNPs.
launcher_creator.py -n 06.mpileup -t 06:00:00 -q normal -w 16 -e -j 06.mpileup_ex.sh
