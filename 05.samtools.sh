#!/usr/bin/env bash 
# Samtools Mapping QC Script Creator
# Chad Smith
#
# Searches for BWA alignments and filters using samtools:
# -q20	Mapping quality score must be >20
# -f 2	Retains alignments in a proper pair
#
# Outputs filtered *.bam files to $BWA_DIR/q20 and $BWA_DIR/q20/pp 
#
# launcher_creator.py must be in $PATH. This creates a script for TACC to launch the jobs. 
#
# BASE          working directory 
# BWA_DIR       dir with *bam alignments
#
# Usage:
# ./05.samtools.sh

BASE=~/scripts/ddrad_pipe/test
BWA_DIR=${BASE}/bwa

if [ ! -d ${BWA_DIR}/q20/pp ]; then mkdir -p ${BWA_DIR}/q20/pp;fi

# Collate names of fq files
find ${BWA_DIR} -maxdepth 1 -iname "*.bam" > bam_raw.txt
nf=`wc -l bam_raw.txt | cut -f1 -d " "`
echo "$nf samples to process."

if [ -f 05.samtools_ex.sh ];then rm 05.samtools_ex.sh;fi
while read line; do
id=`basename ${line} | sed 's;.bam;;'`
echo "samtools view -b -q20 ${line} > ${BWA_DIR}/q20/${id}.bam;\
samtools view -b -f 2 ${BWA_DIR}/q20/${id}.bam > ${BWA_DIR}/q20/pp/${id}.bam" >> 05.samtools_ex.sh
done < bam_raw.txt

launcher_creator.py -n 05.samtools -t 02:00:00 -w 24 -q normal -j 05.samtools_ex.sh
