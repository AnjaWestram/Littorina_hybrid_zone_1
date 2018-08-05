#!/bin/bash -ve

# get allelic depths for each SNP, to be used in cline analysis

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


##############################
#### OUTPUT ALLELIC DEPTH ####
/home/bo1awx/programs/vcftools_0.1.13/bin/vcftools --vcf ANG13c.filt.$taskid.recode.vcf --maf 0.1 --extract-FORMAT-info AD --out CLI01.filt.$taskid

# make column names by repeating each ind twice (once for each allele)
indnames=`head -n 1 CLI01.filt.$taskid.AD.FORMAT | sed "s/CHROM\tPOS\t//g"`
newheader=`echo CHROM POS; for item in $indnames;do echo $item\.1; echo $item\.2; done`
echo $newheader > CLI01.filt.$taskid.alleles1

# data rows
grep -v "CHROM" CLI01.filt.$taskid.AD.FORMAT > CLI01.filt.$taskid.alleles2
sed -i "s/,/\t/g" CLI01.filt.$taskid.alleles2

# combine column names and data rows
cat CLI01.filt.$taskid.alleles1 CLI01.filt.$taskid.alleles2 > CLI01.filt.$taskid.alleles

# rm intermediate files
rm CLI01.filt.$taskid.alleles1
rm CLI01.filt.$taskid.alleles2
