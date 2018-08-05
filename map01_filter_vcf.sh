#!/bin/bash -ve

# filter vcf files for genetic map

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=5:00:00
# current environment settings are used for the job???
#$ -V
#$ -P molecol
#$ -q molecol.q
#$ -t 1000-1987
#$ -tc 20

taskid=${SGE_TASK_ID}

# remove ANG individuals
/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --keep individuals_MAP.txt --vcf ANG_MAP.raw.vcf$taskid --recode --out MAP01.raw.$taskid


#################################
#### GENOTYPE OUTPUT FOR MAP ####
# genotype filters:
/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --vcf MAP01.raw.$taskid.recode.vcf --minGQ 25 --minDP 15 --recode --out MAP01.filt1.$taskid

# variant filters
/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --vcf MAP01.filt1.$taskid.recode.vcf --remove-indels --maf 0.05 --min-alleles 2 --max-alleles 2 --max-missing-count 150 --minQ 20 --recode --out MAP01.filt2.$taskid