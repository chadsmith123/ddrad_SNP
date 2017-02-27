#!/usr/bin/env bash 
# Make Radtag BWA Compliant
# Chad Smith
#
# BWA requires read 1 and read2 have the same header. Replaces _2 in read 2 fastq.gz header 
# to _1 to match read 1. If the fastq header in read2  
# is different from read1 BWA will not work.
# 
# Creates a directory bwa in $LIBS for the BWA compliant reads, and creates soft links in this 
# directory for read 1.
#
# BASE		working directory
# LIB		names of the directories in $BASE containing radtags libraries
#
# Usage:
# ./03.make_bwa_compliant

BASE=~/scripts/ddrad_pipe/test
LIB=(L1 L2)

for i in ${LIB[@]}; do 
 cd ${BASE}/${i}
 if [ ! -d bwa ]; then mkdir bwa; fi

# Copy read 2 to bwa/ and replace last '2' in @fastq header with '1'
 for R2 in `find . -maxdepth 1 -iname "*[ACTG].2.fq.gz"`;do 
   r2_file=`basename $R2`  
   gunzip -c $R2 | sed '/@/ s/2$/1/' | gzip -c > bwa/${r2_file}
 done

# Make soft links for read 1 in bwa/
 for R1 in `find . -maxdepth 1 -iname "*[ACTG].1.fq.gz"`;do 
   r1_file=`basename $R1`  
   r1_dir=`dirname $R1`
   ln -fs ../$r1_file ${r1_dir}/bwa/
 done
done
