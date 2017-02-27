#!/usr/bin/env bash 
# Find and move files
# Chad Smith
#
# Searches filenames in directory $1 for pattern $2 and replaces $3 with $4
if [ -z "$4" ]; then
              echo usage: $0 '<path> <"search string"> <old string> <newstring>'
              exit
          fi
for i in `find $1 -type f -iname "*\$2*"`; do
	originalName="$i"
	# Check that original name is in file search output
	test=`echo -e $originalName | grep $3`				
	if [ -n "$test" ]; then 
       		newName=$(echo -n "$originalName" | sed 's/'$3'/'$4'/')
        	echo "$originalName -> $newName"
        	mv $originalName $newName 
	fi
done
