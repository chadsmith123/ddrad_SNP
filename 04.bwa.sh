#!/usr/bin/env bash
# BWA Script Creator
# Chad Smith
#
# Creates 04.bwa_ex.sh script to map reads on TACC supercomputer.
# The BWA call maps fastq files to $GENOME, discarding secondary alignments (-F 256) and sorting
# the output to *.bam files, one per sample in $BWA_DIR.
#
# launcher_creator.py must be in $PATH. This creates a script for TACC to launch the jobs. 
#
# BASE		working directory
# LIB		names of the directories in $BASE containing radtags libraries
# BWA_DIR	BWA output directory
# GENOME	Path to genome
# TMP		Temporary directory in $SCRATCH for samtools on the cluster
# LIBS		Directory where library(ies) are.
#
# Usage:
# ./04.bwa.sh

BASE=~/scripts/ddrad_pipe/test
BWA_DIR=${BASE}/bwa
GENOME=${WORK}/chad/genomes/atex_v0.1.fa
TMP=${SCRATCH}/tmp
LIBS=${BASE}/L*_radtags

if [ ! -d $BWA_DIR ]; then mkdir $BWA_DIR;fi
if [ ! -d $TMP ]; then mkdir $TMP;fi

# Collate fastq filenames in all libraries
find ${LIBS}/bwa -iname "*.1.fq.gz" > /tmp/all_fq.txt

nfq=`wc -l /tmp/all_fq.txt | cut -f1 -d " "`
count=0
echo "$nfq samples to process in /tmp/all_fq.txt"

if [ -f 04.bwa_ex.sh ];then rm 04.bwa_ex.sh;fi

# Make 04.bwa_ex.sh script
while read line; do
 dir_id=`dirname $line`
 id=`basename $line | sed "s;".1.fq.gz";;"`
 echo "bwa mem $GENOME\
-R '@RG\tID:${id}\tPL:Illumina\tLB:${id}\tSM:${id}'\
${dir_id}/${id}.1.fq.gz ${dir_id}/${id}.2.fq.gz | \
samtools view -bu -F 256 - | \
samtools sort - -T ${TMP}/${id} -o ${BWA_DIR}/${id}.bam" >> 04.bwa_ex.sh 	 

done < /tmp/all_fq.txt

# Create script for running on TACC
launcher_creator.py -n 04.bwa -t 02:00:00 -q development -w 24 -j 04.bwa_ex.sh
