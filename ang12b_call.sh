#!/bin/bash -ve

# variant calling

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=20:00:00
# current environment settings are used for the job???
#$ -V
#$ -P molecol
#$ -q molecol.q
#$ -t 1000-1987
#$ -tc 20

taskid=${SGE_TASK_ID}

/home/bo1awx/programs/samtools-1.3.1/samtools mpileup -b bams.txt -q 20 -Q 20 -L 10000 -t AD,DP -l allcontigs_sorted_common.bed$taskid -u -f Littorina_scaffolded_PacBio_run2_7_Oct_2016_hard_masked.fasta | /home/bo1awx/programs/bcftools-1.3.1/bcftools call -mv -f GQ > ANG_MAP.raw.vcf$taskid

# 1. calculate likelihoods and store as bcf; 2. apply prior and call
# -L: Skip INDEL calling if the average per-input-file depth is above INT
# -u: uncompressed bcf / vcf output
# -t AD: output allelic depth

# -m, --multiallelic-caller:
# alternative model for multiallelic and rare-variant calling designed to overcome known limitations in -c calling model (conflicts with -c) 
# -v, --variants-only:
# output variant sites only 

