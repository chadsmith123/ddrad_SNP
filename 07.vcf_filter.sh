#!/usr/bin/env bash 
# SNP Filter
# Chad Smith
#
# Merges VCF files from mpileup 
# 
# Input are the unmerged VCF files and the output is script 07.vcfmerge_ex.sh to merged VCF files
# on the TACC cluster
#
# MP_OUT	dir of the vcf file output of mpileup
# VCF_UNMERGED	name of the unmerged VCF files
# VCF_PRE	prefix to give the merged VCF file
#
# Requires bcftools
# Usage:
# ./07.vcf_merge.sh

MP_OUT=
VCF_UNMERGED="*raw[0-9].vcf.gz"
VCF_PRE="atex"

# Index VCF files and generate a script to merge them.
cd MP_OUT
ls $VCF_UNMERGED |xargs -I {} bcftools index {}
echo "ls $VCF_UNMERGED |xargs bcftools merge -O z -o ${VCF_PRE}_raw.vcf.gz" > 07.vcfmerge_ex.sh
launcher_creator.py -n 07.SNPvcfmerge -t 02:00:00 -q development -w 24 -e -j 07.vcfmerge_ex.sh
