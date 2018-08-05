#!/bin/bash -ve

# some general filtering

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=01:59:00
# current environment settings are used for the job???
#$ -V
##$ -P molecol
##$ -q molecol.q
#$ -t 1000-1987
##$ -tc 20

taskid=${SGE_TASK_ID}

# remove MAP individuals
/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --keep individuals_ANG.txt --vcf ANG_MAP.raw.vcf$taskid --recode --out ANG13.raw.$taskid

# variant filtering
/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --vcf ANG13.raw.$taskid.recode.vcf --maf 0.01 --minQ 20 --remove-indels --min-alleles 2 --max-alleles 2 --max-missing-count 150 --recode --out ANG13a.filt.$taskid
