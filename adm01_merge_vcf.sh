#!/bin/bash -ve

# merge files from different sets of loci to prepare file for admixture analysis

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=2:00:00
# current environment settings are used for the job???
#$ -V
##$ -P molecol
##$ -q molecol.q

export PERL5LIB=/home/bo1awx/programs/vcftools_0.1.13/perl/

# merge all vcf files
/home/bo1awx/programs/vcftools_0.1.13/bin/vcf-concat ANG14.filt.*.recode.vcf > ADM01.vcf
