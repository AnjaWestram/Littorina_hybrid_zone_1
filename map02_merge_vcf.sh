#!/bin/bash -ve

# merge vcf files for all SNPs and get genotype posterior files fo genetic map

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=8:00:00
# current environment settings are used for the job???
#$ -V
##$ -P molecol
##$ -q molecol.q

export PERL5LIB=/home/bo1awx/programs/vcftools_0.1.13/perl/

# merge all vcf files
/home/bo1awx/programs/vcftools_0.1.13/bin/vcf-concat MAP01.filt2.*.recode.vcf > MAP02.filt2.vcf

# get 012 format
/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --vcf MAP02.filt2.vcf --012 --out MAP02.filt2

# get genotype posteriors format
awk -f map02_vcf2posterior.awk MAP02.filt2.vcf > MAP02.filt2.GP
