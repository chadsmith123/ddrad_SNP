#!/usr/bin/env bash
# Process Radtags
# QC of paired-end Illumina NGS ddRAD data using process_radtags script from 
# STACKS (http://catchenlab.life.illinois.edu/stacks)
#
# Chad Smith
#
# BASE		working directory
# BARCODE	txt file with inline barcodes, one per line
# SEQDIR	location of NGS sequences in fastq.gz format
# OUT_DIR	destination directory of QC'd radtags
# ENZ_1		enzyme 1 used in library prep
# ENZ_2		enzyme 2 used in library prep

BASE=
BARCODE=${BASE}/barcodes_48.txt
SEQDIR=
OUT_DIR=${BASE}/radtags
ENZ_1="nlaIII"
ENZ_2="ecoRI"

if [ ! -d $OUT_DIR ]; then mkdir $OUT_DIR;fi

# -P = paired end reads
# -r rescue barcodes and radrags
# -D output discarded reads
process_radtags -i 'gzfastq' -p $SEQDIR -P -b $BARCODE --inline_null --renz_1 $ENZ_1 --renz_2 $ENZ_2 -r -D -o $OUT_DIR -y fastq
