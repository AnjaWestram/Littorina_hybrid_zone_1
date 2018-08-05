#!/bin/bash -ve

# merge all results of cline analysis

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=01:00:00
# current environment settings are used for the job???
#$ -V

rm -f CLI03_clines.txt

head -n 1 CLI02_cline_snps_ANG-1000.txt > CLI03_clines.txt
cat CLI02_cline_snps_ANG-1*.txt | grep -v Index >> CLI03_clines.txt
