#!/bin/bash -ve

# keep only positions close to mapped ones

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=07:59:00
# current environment settings are used for the job???
#$ -V
##$ -P molecol
##$ -q molecol.q
#$ -t 1000-1987
##$ -tc 20

taskid=${SGE_TASK_ID}

/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --vcf ANG13a.filt.$taskid.recode.vcf --positions ANG13b_positions_to_keep.txt --recode --out ANG13c.filt.$taskid
