#!/bin/bash -ve

# index reference genome

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=07:59:00


/home/bo1awx/programs/bwa-0.7.12/bwa index Littorina_scaffolded_PacBio_run2_7_Oct_2016_hard_masked.fasta
