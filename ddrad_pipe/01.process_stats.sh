#!/bin/bash
BASE=/scratch/02716/cs762/popgen
echo "lib,p_lowquality,p_ambig_tag,p_ambig_bar,retained" > process_radtags_stats.csv
for i in `find ${BASE}/L*_radtags -iname "process_radtags.log"`;do 
	 total=`sed -n '9,13p' $i | grep "Total Sequences" |cut -f2 -d $'\t'`
	 lowq=`sed -n '9,13p' $i | grep "Low Quality" |cut -f2 -d $'\t'`
	 ambig_tag=`sed -n '9,13p' $i | grep "Ambiguous RAD-Tag" |cut -f2 -d $'\t'`
	 retained=`sed -n '9,13p' $i | grep "Retained Reads" |cut -f2 -d $'\t'`
	 ambig_bar=`sed -n '9,13p' $i | grep "Ambiguous Barcodes" |cut -f2 -d $'\t'`
	 p_lowq=`echo "$lowq/$total"|bc`
	 p_ambig_tag=`echo "$ambig_tag/$total"|bc`
	 p_ambig_bar=`echo "$ambig_bar/$total"|bc`
	 p_retained=`echo "$retained/$total"|bc`
	 lib=`dirname $i`
	 
	 echo "$lib,$p_lowq,$p_ambig_tag,$p_ambig_bar,$p_retained,$retained" >>process_radtags_stats.csv
	 done
