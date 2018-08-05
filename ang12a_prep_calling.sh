#!/bin/bash -ve

# prepare files for SNP calling step

# request memory for job (default 6G, max 72G)
#$ -l mem=2G
#$ -l rmem=2G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=02:00:00
# current environment settings are used for the job???
#$ -V

# make file with all bams for calling
rm -f bams.txt
files=`ls *.bam` # all read files
for file1 in $files;do
echo $file1 >> bams.txt
done

# split bed file in order to parallelise calling step
split -d -l 57 -a 3 allcontigs_sorted_common.bed allcontigs_sorted_common.bed1
