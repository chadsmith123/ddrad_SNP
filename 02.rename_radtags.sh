#!/usr/bin/env bash 
# Rename Radtags
# Chad Smith
#
# Renames $prefix to the name of the sample listed in $id using barcodes imbedded in the radtag
# file name.
# !!! Make sure the order of sample names in $id match exactly with the order of barcodes
# in $barcode or radtags will be renamed with the wrong id!
#
# findmove.sh in ddrad_pipe/scripts that must be in your $PATH 
#
# id file example:
# LABienville2_p1_6a
# LABienville2_p1_6b
# LABienville2_p1_6c
# TXDelRio2_p1_6d
#
# corresponding barcode file
# GCATG
# AACCA
# CGATC
# TCGAT
#
# radtag file names with $prefix='sample'
# sample_GCATG.fastq -> LABienville2_p1_6a_GCATG.fastq
# sample_AACCA.fastq -> LABienville2_p1_6b_AACCA.fastq
# sample_CGATC.fastq -> LABienville2_p1_6c_CGATC.fastq
# sample_TCGAT.fastq -> TXDelRio2_p1_6d_TCGAT.fastq
#
# BASE		working directory
# LIBNAME	array of names of the directories in $BASE containing radtags libraries
# prefix	common prefix of radtags files
# id		id file path. Put one in each $LIBNAME.
# barcodes	barcodes file path
#
# Usage:
# ./02.rename_radtags.sh

BASE=~/scripts/ddrad_pipe/test
LIBNAME=(L1 L2 L3 L4 L5 L6 L7 L8 L9)

prefix="sample"
barcodes=barcodes.txt
id=id.txt

function findmove_radtags {
 # Sanity check that $id and $barcodes exist and have the same number of lines.
 if [ ! -f "$id" ] | [ ! -f "$barcodes" ]; then echo "$id or $barcodes files not found. Exiting."; exit; fi
 
 wc_bar=`wc -l $barcodes |awk '{print $1}'`
 wc_id=`wc -l $id | awk '{print $1}'`
 
 if [ "$wc_bar" != "$wc_id" ]; then echo "Barcode file has $wc_bar entries and id file in `pwd`\
 has $wc_id entries. Script stopped."; exit
  else
 
 # Cleanup previous runs of script
 if [ -f /tmp/sample.txt ]; then rm /tmp/sample.txt; fi
 if [ -f rename_samples.log ]; then rm rename_samples.log; fi
 
 # Make rename_samples.sh script from $id and $barcodes and execute
 for i in `cat $barcodes`; do echo $prefix >> /tmp/sample.txt; done
 
 sed -e 's/^/"/' -e 's/$/"/' $barcodes > /tmp/bar2.txt
 paste -d " " /tmp/bar2.txt /tmp/sample.txt $id > /tmp/paste.txt
 sed -e 's;^;findmove.sh . ;' -e 's;$; >> rename_samples.log;' /tmp/paste.txt > rename_samples.sh
 chmod +x rename_samples.sh
 ./rename_samples.sh
 fi
}

# Run
for i in ${LIBNAME[@]}; do 
 cd ${BASE}/${i}
 findmove_radtags
done
