# ddrad_SNP
## Pipeline for Double-Digest Restriction Associated DNA Sequencing

This pipeline processes paired-end Illumina reads from a ddRAD library prep and outputs SNPs using samtools mpileup. The steps are:

1) Run quality control of radtags 
2) Rename generic file names (eg sample*.fq.gz) with a user-provided name in the file "id.txt"
3) Change the header of fastq sequences to be compliant with BWA
4) Generate the 04.bwa_ex.sh script to align reads on the TACC cluster 
5) Filter alignments 
6) Generatea the script 06.mpileup_ex.sh script to call SNPS on the TACC cluster 
7) Generate the 07.vcfmerge_ex.sh script to merge vcf files on the TACC cluster 
8) Filter SNPs for population genomic analysis

## Dependencies
- STACKS
- bc
- bwa mem
- samtools
- bcftools
- vcftools 

You must add the scripts/ directory to your PATH. In your bash config file (usually ~/.bashrc or ~/.bash_profile)
EXPORT PATH=${PATH}:<name of directory with ddrad_SNP/scripts>

## Usage
1) Create a directory to put the output from the pipeline. This is referred to as $BASE in the scripts
2) Create a file barcodes.txt with a list of the internal barcodes.
3) Make a directory for the radtags in $BASE and put a text file id.txt with the IDs of your samples. These MUST be ordered to correspond with the barcodes in the barcodes.txt file.
For example, if the first barcode in barcode.txt is "ACTGC", the the first sample name in id.txt must be the sample that was tagged with the "ACTGC" barcode.
4) Edit the scripts to change the parameters to your liking.
5) cd to the $BASE directory and run the scripts in order.

At the end of the pipeline $BASE will have:

- a dir with the radtags 
- a dir with the alignments
- a dir with mpileup output
