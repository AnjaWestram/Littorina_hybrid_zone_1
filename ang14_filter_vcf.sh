#!/bin/bash -ve

# filtering used before admixture analysis, LD analysis, BayeScan

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

# filter
/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --vcf ANG13c.filt.$taskid.recode.vcf --minGQ 20 --minDP 15 --recode --out ANG14.filt.temp.$taskid

/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --vcf ANG14.filt.temp.$taskid.recode.vcf --maf 0.01 --minQ 20 --max-missing-count 150 --recode --out ANG14.filt.$taskid


##########################
#### OUTPUT GENOTYPES ####
/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --vcf ANG14.filt.$taskid.recode.vcf --extract-FORMAT-info GT --out ANG14.filt.$taskid

