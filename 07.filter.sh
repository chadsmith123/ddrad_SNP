#!/usr/bin/env bash 
# SNP Filter
# Chad Smith
#
# Merges VCF files from mpileup and performs QC using vcftools.
# 
# The order of operation is: 
# 1) Generates a script to merge VCFs on the cluster
# 2) Filter by sequencing depth and genotyping quality
# 3) Remove samples (duplicates, failed samples, etc)
# 4) Thin sites to 1/100k KB
# 5) Generate .012 and X/Y formatted files for input to adegenet
#
# Input are the unmerged VCF files and the output is a filtered VCF file, a VCF with a subset of 
# 1000k sites, and 012 formatted files for adegenet
#
# VCF_UNMERGED	name of the unmerged VCF files
# VCF_PRE	prefix to give the merged VCF file
# RM		text file with a list of any samples to remove from VCF file
# minGQ		minimum genotyping quality to retain a site in a given sample
# minDP		minimum depth of sequencing to retain a site in a given sample
# pMISS		minimum % of samples with a site to retain it. 
#		0=allow sites completely missing, 1=no missing data allowed
#
# Requires vcftools, bcftools
# Usage:
# ./07.filter.sh

VCF_UNMERGED="*raw[0-9].vcf.gz"
VCF_PRE="atex"
RM=rm_samples.txt
minGQ=20
minDP=15
pMISS=0.8

# Index VCF files and generate a script to merge them.
ls $VCF_UNMERGED |xargs -I {} bcftools index {}
echo "ls $VCF_UNMERGED |xargs bcftools merge -O z -o ${VCF_PRE}_raw.vcf.gz" > 07.filter_merge_ex.sh
launcher_creator.py -n 07.SNPfilter_merge -t 02:00:00 -q development -w 24 -e -j 07.filter_merge_ex.sh

vcftools --gzvcf ${VCF_PRE}_raw.vcf.gz --max-missing ${pMISS} --minDP ${minDP} --minGQ ${minGQ} --recode --recode-INFO-all --stdout > ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}.vcf

# Remove samples (duplicates, failed library preps, etc)
vcftools --vcf ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}.vcf --remove $RM --recode --recode-INFO-all --stdout > ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}_rm.vcf

# Thin SNPs so there are no SNPs w/in 100KB. This reduces linkage between SNPs.
vcftools --vcf ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}_rm.vcf --thin 100000 --recode --recode-INFO-all --stdout > ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}_rm_100kb.vcf

# Take a subset of 1000 SNPs
vcftools --vcf ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}_rm_100kb.vcf --site-mean-depth --stdout| cut -f1,2 | shuf -n 10 > ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}_rm_100kb.1k
bcftools view -T ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}_rm_100kb.1k ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}_rm_100kb.vcf > ${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}_rm_100kb_1k.vcf

## Generate genotypes file compatible with adegenet
# Create 012 allele coding file
PREF=${VCF_PRE}_gq${minGQ}_DP${minDP}_p${pMISS}_rm_100kb_1k
vcftools --vcf ${PREF}.vcf --012 --out ${PREF} 
cut -f2- ${PREF}.012 > /tmp/${PREF}.012.tmp
paste ${PREF}.012.indv /tmp/${PREF}.012.tmp > ${PREF}.012.csv
sed -i 's;\t;_;' ${PREF}.012.pos

# Create X/X format from 012 coding file
cut -f 2- ${PREF}.012 |sed -e 's;-1;NA;g' -e 's;0;0/0;g' -e 's;1;0/1;g' -e 's;2;1/1;g' > /tmp/${PREF}_slash
paste ${PREF}.012.indv /tmp/${PREF}_slash > ${PREF}_slash.csv
