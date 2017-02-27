#!/usr/bin/env bash
# Process Radtags
# Chad Smith
#
# QC of paired-end Illumina NGS ddRAD data using process_radtags script from 
# STACKS (http://catchenlab.life.illinois.edu/stacks)
# Output is gzip'd fastq files to $OUT_DIR and input is raw Illumina fastq files from $SEQ_DIR.
# 
# SEQDIR	location of NGS sequences in fastq.gz format
# BASE		working directory
# LIBNAME	name to give the radtags library. A directory in $BASE will be made with this name.
# BARCODE	txt file with inline barcodes, one per line
# ENZ_1		enzyme 1 used in library prep
# ENZ_2		enzyme 2 used in library prep
#
# Usage:
# ./01.process_radtags.sh

BASE=
SEQDIR=
LIBNAME=
BARCODE=${BASE}/barcodes_48.txt
ENZ_1="nlaIII"
ENZ_2="ecoRI"

OUT_DIR=${BASE}/{$LIBNAME}
function process_radtags_log () {
 # Process_radtag Statistics Generator
 # Chad Smith
 #
 # Parses the $LOG from Stacks process_radtags and outputs summary statistics to the file 
 # process_radtags_stats.csv 
 #
 # lib           name of the ddRAD library
 # LOG		 name of the radtag log file 
 # lowq          number of low quality sequences 
 # ambig_tag     number of ambiguous radtags
 # retained      number of reads retained after filtering
 # ambig_bar     number of ambigous barcodes
 # p_X           as above but converted into a proportion
 #
 # Requires bc to do the math: https://www.gnu.org/software/bc/manual/html_mono/bc.html
 #
 # LOG		 path to log file
 #
 # Usage
 # process_radtags <lib_name>

 LOG='process_radtags.log'

 lib=$1

 total=`grep "^Total Sequences" $LOG |cut -f2 -d $'\t'`
 lowq=`grep "^Low Quality" $LOG |cut -f2 -d $'\t'`
 ambig_tag=`grep "^Ambiguous RAD-Tag" $LOG |cut -f2 -d $'\t'`
 retained=`grep "^Retained Reads" $LOG |cut -f2 -d $'\t'`
 ambig_bar=`grep "^Ambiguous Barcodes" $LOG |cut -f2 -d $'\t'`

 p_lowq=`echo "$lowq/$total"|bc`
 p_ambig_tag=`echo "$ambig_tag/$total"|bc`
 p_ambig_bar=`echo "$ambig_bar/$total"|bc`
 p_retained=`echo "$retained/$total"|bc`
         
 echo "lib,p_lowquality,p_ambig_tag,p_ambig_bar,retained" > process_radtags_stats.csv
 echo "$lib,$p_lowq,$p_ambig_tag,$p_ambig_bar,$p_retained,$retained" >>process_radtags_stats.csv
}

# Run
if [ ! -d $OUT_DIR ]; then mkdir $OUT_DIR;fi

# -P = paired end reads
# -r rescue barcodes and radrags
# -D output discarded reads
process_radtags -i 'gzfastq' -p $SEQDIR -P -b $BARCODE --inline_null --renz_1 $ENZ_1 --renz_2 $ENZ_2 -r -D -o $OUT_DIR -y gzfastq
cd $OUT_DIR; process_radtags_log $LIBNAME
