#!/bin/bash -ve

# run LD / dispersal analysis

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=167:59:00
#$ -V
##$ -P molecol
##$ -q molecol.q
##$ -t 1000-1987
##$ -tc 20

taskid=${SGE_TASK_ID}

module add apps/R/3.3.1

Rscript --vanilla ang15_LD.R