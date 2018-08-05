#!/bin/bash -ve

# get all positions of biallelic SNPs (to make a file from which only those close to mapped positions will be selected)

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=07:59:00
# current environment settings are used for the job???
#$ -V
#$ -P molecol
#$ -q molecol.q

for file1 in ANG13a.filt.1*.recode.vcf;do
grep -v "#" $file1 | cut -f1,2 >> ANG13b_all_snps.txt
done

sort ANG13b_all_snps.txt | uniq > ANG13b_all_snps_unique.txt

